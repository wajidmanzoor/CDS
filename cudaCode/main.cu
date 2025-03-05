#include "./inc/graph.h"
#include "./inc/motif.h"

#include "./src/gpuMemoryAllocation.cc"
#include "./src/helpers.cc"


#define CUDA_CHECK_ERROR(kernelName) { \
    cudaError_t err = cudaGetLastError(); \
    if (err != cudaSuccess) { \
        printf("CUDA Error in kernel %s, file %s at line %d: %s\n", kernelName, __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
    } \
}

bool fileExists(const string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

void writeOrAppend(const string& filename, const string& data) {
    ofstream file;
    
    // Check if the file exists
    if (fileExists(filename)) {
        // Open the file in append mode if it exists
        file.open(filename, ios::app);
    } else {
        // Open the file in write mode if it doesn't exist
        file.open(filename);
    }
    
    if (file.is_open()) {
        file << data << endl;
        file.close();
    } else {
        cerr << "Unable to open the file." << endl;
    }
}

int main(int argc, const char * argv[]) {
    if (argc != 6) {
        cout << "Server wrong input parameters!" << endl;
        exit(1);
    }

    string filepath = argv[1]; // Path to the graph file. The graph should be represented as an adjacency list with space separators
    string motifPath = argv[2]; //Path to motif file. The motif should be represented as edge list with space sperators.
    ui k = atoi(argv[3]);
    ui pSize = atoi(argv[4]);
    ui cpSize = atoi(argv[5]);

    Graph G = Graph(filepath);

    Motif M = Motif(motifPath);

    vector<ui> *listingOrder;
    
    G.getListingOrder(&listingOrder);
    memoryAllocationGraph(deviceGraph,G);
    memoryAllocationDAG(deviceDAG,G.n;G.m);

    // THIS PART IS RELATED TO GENERATING DAG

    ui *listOrder;
    chkerr(cudaMalloc((void**)&(listOrder), n * sizeof(ui)));
    chkerr(cudaMemcpy(listOrder, listingOrder->data(), n * sizeof(ui), cudaMemcpyHostToDevice));

    // Get out degree in DAG
    generateDegreeDAG<<<BLK_NUMS, BLK_DIM>>>(deviceGraph,deviceDAG,listOrder,G.n, G.m,TOTAL_WARPS);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Generate Degree of DAG");

    //copy out degree to offset
    chkerr(cudaMemset(deviceDAG.neighbor, 0, sizeof(ui)));
    chkerr(cudaMemcpy(deviceDAG.offset + 1,deviceDAG.degree , (G.n) * sizeof(ui), cudaMemcpyDevicetoDevice));

    // cummulative sum offset
    thrust::inclusive_scan(thrust::device_ptr < ui > (deviceDAG.offset), thrust::device_ptr < ui > (deviceDAG.offset + G.n + 1), thrust::device_ptr < ui > (deviceDAG.offset ));

    // Write neighbors of DAG
    generateNeighborDAG<<<BLK_NUMS, BLK_DIM>>>(deviceGraph,deviceDAG,listOrder,G.n, G.m,TOTAL_WARPS);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Generate Neighbors of DAG");
    chkerr(cudaFree(listOrder));


    // THIS PART IS ABOUT CLIQUE LISTING ALGORITHM

    int maxDegree =0;
    ui maxBitMask = memoryAllocationlevelData(levelData,k,pSize,cpSize,maxDegree,TOTAL_WARPS);
    int level =0;
    int iterK = k;

    ui *labels;
    chkerr(cudaMalloc((void**)&(labels), n * sizeof(ui)));

    thrust::device_ptr<ui> dev_labels(labels);
    thrust::fill(dev_labels, dev_labels + n, iterK);

    chkerr(cudaMemcpy(G.degree, graph.degree->data(), n * sizeof(ui), cudaMemcpyHostToDevice));

    //TODO SHARED MEMORY 
    listIntialCliques<<<BLK_NUMS, BLK_DIM>>>(deviceDAG, levelData, label,iterK, G.n, G.m,psize,cpSize,maxBitMask,level,TOTAL_WARPS);
    CUDA_CHECK_ERROR("Generate Intial Partial Cliques");

    iterK --;
    level ++;
    ui offsetPartitionSize = ((psize/(k-1)) + 1);
    thrust::transform(thrust::device, thrust::make_counting_iterator(0), thrust::make_counting_iterator(TOTAL_WARPS), leveldata.temp + 1,
                      [leveldata.offsetPartition, leveldata.count,offsetPartitionSize] __device__ (int i) {
                          int task_count = leveldata.count[i];
                          return (task_count > 0) ? leveldata.offsetPartition[i * offsetPartitionSize + task_count] : 0;
                      });

    thrust::inclusive_scan(thrust::device, leveldata.temp, leveldata.temp + TOTAL_WARPS + 1, leveldata.temp);
    thrust::inclusive_scan(thrust::device, leveldata.count, leveldata.count + TOTAL_WARPS + 1, leveldata.count);

    flushParitions<<<BLK_NUMS, BLK_DIM>>>( deviceDAG, levelData,psize,cpSize,maxBitMask, level,TOTAL_WARPS);
    CUDA_CHECK_ERROR("Flush Partition data structure");

    int totalTasks;
    chkerr(cudaMemcpy(&totalTasks,leveldata.count + TOTAL_WARPS, sizeof(ui), cudaMemcpyDeviceToHost));

    
    while(iterK > 2 ){

        //TODO: need to think a solution for label . Maybe use it masks will see 
        thrust::fill(dev_labels, dev_labels + n, iterK);
        listMidCliques(deviceDAG, levelData,label,k,iterK, G.n, G.m,psize,cpSize,maxBitMask,totalTasks,level,TOTAL_WARPS);
        CUDA_CHECK_ERROR("Generate Mid Partial Cliques");

        thrust::transform(thrust::device, thrust::make_counting_iterator(0), thrust::make_counting_iterator(TOTAL_WARPS), leveldata.temp + 1,
                      [leveldata.offsetPartition, leveldata.count,offsetPartitionSize] __device__ (int i) {
                          int task_count = leveldata.count[i];
                          return (task_count > 0) ? leveldata.offsetPartition[i * offsetPartitionSize + task_count] : 0;
                      });

        thrust::inclusive_scan(thrust::device, leveldata.temp, leveldata.temp + TOTAL_WARPS + 1, leveldata.temp);
        thrust::inclusive_scan(thrust::device, leveldata.count, leveldata.count + TOTAL_WARPS + 1, leveldata.count);

        flushParitions<<<BLK_NUMS, BLK_DIM>>>( deviceDAG, levelData,psize,cpSize,k,maxBitMask,level,TOTAL_WARPS);
        CUDA_CHECK_ERROR("Flush Partition data structure");

        chkerr(cudaMemcpy(&totalTasks,leveldata.count + TOTAL_WARPS, sizeof(ui), cudaMemcpyDeviceToHost));
        iterK --; 
        level ++;
    }
    ui t; 
    //TODO: decide the total number cliques
    chkerr(cudaFree(label));
    freeLevelPartitionData(levelData);
    
    memoryAllocationTrie(cliqueData, t, k);
    int totalCliques;

    chkerr(cudaMalloc((void**)&totalCliques,sizeof(ui)));
    chkerr(cudaMemset(totalCliques,0,sizeof(ui)));

    if(iterK == 2){
        writeFinalCliques<<<BLK_NUMS, BLK_DIM>>>(deviceGraph, deviceDAG,  levelData,cliqueData,totalCliques,k,iterK, n, m,psize,cpSize, maxBitMask,totalTasks, level, TOTAL_WARPS);
        CUDA_CHECK_ERROR("Generate Full Cliques");
        
    }
    freeLevelData(levelData);
    freeDAG(deviceDAG);

    //TODO:  CLIQUE CORE DECOMPOSE

    //TODO: LOCATE CORE

    //TODO: LISTING AGAIN

    //TODO: EDGE PRUNING

    //TODO: COMPONENT DECOMPOSE

    //TODO: DYNAMIC CORE EXACT






     
    
    
    freeGraph(deviceGraph);
    delete G;
    delete M;
    return 0;
}
