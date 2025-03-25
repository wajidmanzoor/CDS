#include "./inc/common.h"
#include "./inc/graph.h"

#include "./utils/cuda_utils.cuh"
#include "./inc/gpuMemoryAllocation.cuh"
#include "./inc/helpers.cuh"


void generateDAG(const Graph& graph,deviceGraphPointers& deviceGraph,deviceDAGpointer& deviceDAG , vector<ui> listingOrder){
    
    memoryAllocationDAG(deviceDAG, graph.n, graph.m);

    ui *listOrder;
    chkerr(cudaMalloc((void**)&(listOrder), graph.n * sizeof(ui)));
    chkerr(cudaMemcpy(listOrder, listingOrder.data(), graph.n * sizeof(ui), cudaMemcpyHostToDevice));

    // Get out degree in DAG
    generateDegreeDAG<<<BLK_NUMS, BLK_DIM>>>(deviceGraph, deviceDAG, listOrder, graph.n, graph.m, TOTAL_WARPS);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Generate Degree of DAG");


    //copy out degree to offset
    chkerr(cudaMemset(deviceDAG.offset, 0, sizeof(ui)));
    chkerr(cudaMemcpy(deviceDAG.offset + 1, deviceDAG.degree, (graph.n) * sizeof(ui), cudaMemcpyDeviceToDevice));

    // cummulative sum offset
    thrust::inclusive_scan(thrust::device_ptr<ui>(deviceDAG.offset), thrust::device_ptr<ui>(deviceDAG.offset + graph.n + 1), thrust::device_ptr<ui>(deviceDAG.offset));


    // Write neighbors of DAG
    size_t sharedMemoryGenDagNeig =  WARPS_EACH_BLK * sizeof(ui);
    generateNeighborDAG<<<BLK_NUMS, BLK_DIM,sharedMemoryGenDagNeig>>>(deviceGraph, deviceDAG, listOrder, graph.n, graph.m, TOTAL_WARPS);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Generate Neighbor of DAG");

    chkerr(cudaFree(listOrder));

}

ui listAllCliques(const Graph& graph,deviceGraphPointers& deviceGraph,deviceDAGpointer& deviceDAG,cliqueLevelDataPointer levelData, ui k, ui pSize,ui  cpSize, ui t){
    thrust::device_ptr<ui> dev_degree(deviceDAG.degree);


    auto max_iter = thrust::max_element(dev_degree, dev_degree + graph.n);
    int maxDegree = *max_iter;

    ui maxBitMask = memoryAllocationlevelData(levelData, k, pSize, cpSize, maxDegree, TOTAL_WARPS);

    int level = 0;
    int iterK = k;

    ui *labels;
    chkerr(cudaMalloc((void**)&(labels), (graph.n * TOTAL_WARPS) * sizeof(ui)));
    thrust::device_ptr<ui> dev_labels(labels);
    thrust::fill(dev_labels, dev_labels + graph.n*TOTAL_WARPS, iterK);

    chkerr(cudaMemcpy(deviceGraph.degree, graph.degree.data(), graph.n * sizeof(ui), cudaMemcpyHostToDevice));
    chkerr(cudaMemset(levelData.partialCliquesPartition, 0,  (TOTAL_WARPS * pSize)* sizeof(ui)));

    size_t sharedMemoryIntialClique =  WARPS_EACH_BLK * sizeof(ui);
    listIntialCliques<<<BLK_NUMS, BLK_DIM,sharedMemoryIntialClique>>>(deviceDAG, levelData, labels, iterK, graph.n, graph.m, pSize, cpSize, maxBitMask, level, TOTAL_WARPS);
    CUDA_CHECK_ERROR("Generate Intial Partial Cliques");
    ui partialSize = TOTAL_WARPS * pSize;
    //ui candidateSize = TOTAL_WARPS * cpSize;
    ui offsetSize = ((pSize / (k - 1)) + 1) * TOTAL_WARPS;

  
    ui offsetPartitionSize = ((pSize / (k-1)) + 1);
    createLevelDataOffset(levelData, offsetPartitionSize, TOTAL_WARPS);

    flushParitions<<<BLK_NUMS, BLK_DIM>>>(deviceDAG, levelData, pSize,cpSize,k, maxBitMask, level, TOTAL_WARPS);
    CUDA_CHECK_ERROR("Flush Partition data structure");

    iterK--;
    level++;

    int totalTasks;
    chkerr(cudaMemcpy(&totalTasks, levelData.count + TOTAL_WARPS, sizeof(ui), cudaMemcpyDeviceToHost));
    size_t sharedMemoryMid =  WARPS_EACH_BLK * sizeof(ui);

    while(iterK > 2) {
        thrust::device_ptr<ui> dev_labels(labels);
        thrust::fill(dev_labels, dev_labels + graph.n*TOTAL_WARPS, iterK);
        chkerr(cudaMemset(levelData.count, 0, (TOTAL_WARPS + 1) * sizeof(ui)));
        chkerr(cudaMemset(levelData.temp, 0, (TOTAL_WARPS + 1) * sizeof(ui)));
        chkerr(cudaMemset(levelData.offsetPartition, 0,  (offsetSize)* sizeof(ui)));
        chkerr(cudaMemset(levelData.validNeighMaskPartition,0, (partialSize * maxBitMask) * sizeof(ui)));

        listMidCliques<<<BLK_NUMS, BLK_DIM,sharedMemoryMid>>>(deviceDAG, levelData, labels, k, iterK, graph.n, graph.m, pSize, cpSize, maxBitMask, totalTasks, level, TOTAL_WARPS);
        CUDA_CHECK_ERROR("Generate Mid Partial Cliques");

        createLevelDataOffset(levelData, offsetPartitionSize, TOTAL_WARPS);

        chkerr(cudaMemset(levelData.offset,0,offsetSize*sizeof(ui)));
        chkerr(cudaMemset(levelData.validNeighMask,0,partialSize*maxBitMask*sizeof(ui)));

        flushParitions<<<BLK_NUMS, BLK_DIM>>>(deviceDAG, levelData, pSize,cpSize,k, maxBitMask, level, TOTAL_WARPS);
        CUDA_CHECK_ERROR("Flush Partition data structure");

        iterK--;
        level++;
    }

    //ui t = 10; // Make it ud
    //TODO: decide the total number cliques and Free Level Data p 1
    chkerr(cudaFree(labels));    
    memoryAllocationTrie(cliqueData, t, k);

    ui *totalCliques;
    chkerr(cudaMalloc((void**)&totalCliques, sizeof(ui)));
    chkerr(cudaMemset(totalCliques, 0, sizeof(ui)));
    size_t sharedMemoryFinal =  WARPS_EACH_BLK * sizeof(ui);


    chkerr(cudaMemset(cliqueData.status, 0, t * sizeof(ui)));
    chkerr(cudaMemset(cliqueData.trie, 0, t * k * sizeof(ui)));
    if(iterK == 2) {
        writeFinalCliques<<<BLK_NUMS, BLK_DIM,sharedMemoryFinal>>>(deviceGraph, deviceDAG, levelData, cliqueData, totalCliques, k, iterK, graph.n, graph.m, pSize, cpSize, maxBitMask, t,totalTasks, level, TOTAL_WARPS);
        CUDA_CHECK_ERROR("Generate Full Cliques");
    }


    freeLevelData(levelData);
    freeLevelPartitionData(levelData);
    freeDAG(deviceDAG);

    ui tt;
    chkerr(cudaMemcpy(&tt, totalCliques, sizeof(ui), cudaMemcpyDeviceToHost));
    cout<<endl<<"total cliques "<<tt<<endl;

    size_t sharedMemorySort =  2*k*WARPS_EACH_BLK * sizeof(ui);
    sortTrieData<<<BLK_NUMS, BLK_DIM,sharedMemorySort>>>(deviceGraph, cliqueData, tt,t, k, TOTAL_THREAD);
    CUDA_CHECK_ERROR("Sort Trie Data Structure");

    cudaFree(totalCliques)
    return tt;
}

void cliqueCoreDecompose(const Graph& graph,deviceGraphPointers& deviceGraph,deviceCliquesPointer& cliqueData, ui &maxCore, double &maxDensity, int &argmax, ui &coreSize, ui &coreTotalCliques, ui glBufferSize, ui k, ui t, ui tt){
    ui level = 0;
    ui count = 0;
    ui *globalCount = NULL;
    ui *bufTails  = NULL;
    ui *glBuffers = NULL;

    chkerr(cudaMalloc((void**)&(globalCount), sizeof(ui)));
    chkerr(cudaMalloc((void**)&(bufTails), BLK_NUMS*sizeof(ui)));
    chkerr(cudaMalloc((void**)&(glBuffers), BLK_NUMS*glBufferSize*sizeof(ui)));
    chkerr(cudaMemset(globalCount, 0, sizeof(ui)));
    cudaDeviceSynchronize();

    chkerr(cudaMalloc((void**)&coreDensity, (graph.n+1)*sizeof(double)));

    ui *removedVerticies, *remainingCliques;
    chkerr(cudaMalloc((void**)&removedVerticies, (graph.n+1)*sizeof(ui)));
    chkerr(cudaMalloc((void**)&remainingCliques, (graph.n+1)*sizeof(ui)));


//    chkerr(cudaMemset(glBuffers, 0, BLK_NUMS*glBufferSize*sizeof(ui)));

    chkerr(cudaMemcpy(deviceGraph.cliqueCore, deviceGraph.cliqueDegree, graph.n * sizeof(ui), cudaMemcpyDeviceToDevice));

    //density of full graph
    thrust::device_vector<ui> dev_vec1(cliqueData.status, cliqueData.status + t);
    ui currentCliques = thrust::reduce(dev_vec1.begin(), dev_vec1.end(), 0, thrust::plus<ui>());
    double d = static_cast<double>(currentCliques) / (graph.n - count);
    chkerr(cudaMemcpy(coreDensity+level, &d, sizeof(double), cudaMemcpyHostToDevice));

    while(count < graph.n){
        cudaMemset(bufTails, 0, sizeof(unsigned int)*BLK_NUMS);

        // Select nodes whoes current degree is level, that means they should be removed as part of the level core 
        selectNodes<<<BLK_NUMS, BLK_DIM>>>(deviceGraph, bufTails, glBuffers, glBufferSize, graph.n, level);
        cudaDeviceSynchronize();
        
        //Total number of verticies in buffer
        thrust::device_vector<ui> dev_vec(bufTails, bufTails + BLK_NUMS);
        ui sum = thrust::reduce(dev_vec.begin(), dev_vec.end(), 0, thrust::plus<ui>());

        //Bases on total vertices device to either use Warp or Block to process one vertex and its cliques
        if(true){
            processNodesByWarp<<<BLK_NUMS, BLK_DIM>>>(deviceGraph, cliqueData , bufTails, glBuffers, globalCount, glBufferSize, graph.n, level, k, tt);
            cudaDeviceSynchronize();
        }else{
            processNodesByBlock<<<BLK_NUMS, BLK_DIM>>>(deviceGraph, cliqueData , bufTails, glBuffers, globalCount, glBufferSize, graph.n, level, k, tt);
            cudaDeviceSynchronize();

        }

        chkerr(cudaMemcpy(&count, globalCount, sizeof(unsigned int), cudaMemcpyDeviceToHost));    
        level++;
        currentCliques = thrust::reduce(dev_vec1.begin(), dev_vec1.end(), 0, thrust::plus<ui>());
        d = static_cast<double>(currentCliques) / (graph.n - count);

        chkerr(cudaMemcpy(coreDensity+level, &d, sizeof(double), cudaMemcpyHostToDevice));
        chkerr(cudaMemcpy(removedVerticies+level , &count, sizeof(ui), cudaMemcpyHostToDevice));
        chkerr(cudaMemcpy(remainingCliques+level , &currentCliques, sizeof(ui), cudaMemcpyHostToDevice));
        cudaDeviceSynchronize();

    }
    graph.kmax = level-1;
    maxCore = level-1;
    thrust::device_ptr<double> coreDensity_ptr(coreDensity);
    thrust::device_ptr<double> max_iter1 = thrust::max_element(coreDensity_ptr, coreDensity_ptr + level);
    maxDensity = *max_iter1;
    argmax = max_iter1 - coreDensity_ptr;

    chkerr(cudaMemcpy(&coreSize, removedVerticies + argmax, sizeof(ui), cudaMemcpyDeviceToHost));
    coreSize = graph.n - coreSize;
    chkerr(cudaMemcpy(&coreTotalCliques, remainingCliques + argmax, sizeof(ui), cudaMemcpyDeviceToHost));

    cudaFree(globalCount);
    cudaFree(bufTails);
    cudaFree(glBuffers);
    cudaFree(coreDensity);
    cudaFree(removedVerticies);
    cudaFree(remainingCliques);

}
int main(int argc, const char * argv[]) {
    if (argc != 7) {
        cout << "Server wrong input parameters!" << endl;
        exit(1);
    }

    string filepath = argv[1]; // Path to the graph file. The graph should be represented as an adjacency list with space separators
    string motifPath = argv[2]; //Path to motif file. The motif should be represented as edge list with space sperators.
    ui k = atoi(argv[3]);
    ui pSize = atoi(argv[4]);
    ui cpSize = atoi(argv[5]);
    ui glBufferSize = atoi(argv[6]);

    // Need better way to do this
    ui t=10;

    Graph graph = Graph(filepath);

    //Motif M = Motif(motifPath);

    vector<ui> listingOrder;
    listingOrder.resize(graph.n);
    graph.getListingOrder(listingOrder);

    memoryAllocationGraph(deviceGraph, graph);

    // GENERATES DAG
    generateDAG(graph, deviceGraph, deviceDAG,listingOrder);

    // CLIQUE LISTING ALGORITHM
    ui tt = listAllCliques(graph, deviceGraph, deviceDAG, levelData, k, pSize, cpSize,t);


    // Debug start 
    int *h_cliques,*status;
    h_cliques = new int[t*k];
    status = new int[t];
    chkerr(cudaMemcpy(h_cliques, cliqueData.trie, k * t * sizeof(ui), cudaMemcpyDeviceToHost));
    chkerr(cudaMemcpy(status, cliqueData.status, t * sizeof(ui), cudaMemcpyDeviceToHost));
    cout<<endl;
    for(int i =0;i<k;i++){
      cout<<endl<<"CL "<<i<<"  ";
      for(int j =0;j<t;j++){
        cout<<h_cliques[i*t+j]<<" ";
      }
    }
    cout<<endl<<"stat  ";
    for(int i = 0; i < t; i++) {
        cout << status[i] << " ";

    }

    ui *h_cdegree;
    h_cdegree = new ui[graph.n];
    chkerr(cudaMemcpy(h_cdegree, deviceGraph.cliqueDegree, graph.n* sizeof(ui), cudaMemcpyDeviceToHost));
    cout<<endl;
    for(int i = 0; i < graph.n; i++) {
        cout <<i<<" ";
    } 
    cout<<endl;
    for(int i = 0; i < graph.n; i++) {
        cout << h_cdegree[i] << " ";
    }

    //Debug end

    //TODO:  CLIQUE CORE DECOMPOSE
    ui coreSize, coreTotalCliques,maxCore;
    double maxDensity;
    int argmax;
    cliqueCoreDecompose(graph,deviceGraph,cliqueData,maxCore, maxDensity, argmax, coreSize, coreTotalCliques,glBufferSize, k,  t, tt);


    

    //TODO: LOCATE CORE

    //TODO: LISTING AGAIN

    //TODO: EDGE PRUNING

    //TODO: COMPONENT DECOMPOSE

    //TODO: DYNAMIC CORE EXACT
    freTrie(cliqueData);
    freeGraph(deviceGraph);
    return 0;
}