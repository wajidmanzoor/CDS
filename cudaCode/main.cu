#include "./inc/common.h"
#include "./inc/graph.h"
#include "./inc/motif.h"

#include "./utils/cuda_utils.cuh"

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

    Graph graph = Graph(filepath);

    //Motif M = Motif(motifPath);

    vector<ui> listingOrder;
    listingOrder.resize(graph.n);
    graph.getListingOrder(listingOrder);

    memoryAllocationGraph(deviceGraph, graph);
    memoryAllocationDAG(deviceDAG, graph.n, graph.m);

    // THIS PART IS RELATED TO GENERATING DAG

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


    // THIS PART IS ABOUT CLIQUE LISTING ALGORITHM

    thrust::device_ptr<ui> dev_degree(deviceDAG.degree);


    auto max_iter = thrust::max_element(dev_degree, dev_degree + graph.n);

    int maxDegree = *max_iter;

    cout<<"max deg "<<maxDegree<<endl;
    ui maxBitMask = memoryAllocationlevelData(levelData, k, pSize, cpSize, maxDegree, TOTAL_WARPS);

    cout<<"max bits "<<maxBitMask<<endl;
    int level = 0;
    int iterK = k;

    ui *labels;
    chkerr(cudaMalloc((void**)&(labels), (graph.n * TOTAL_WARPS) * sizeof(ui)));
    thrust::device_ptr<ui> dev_labels(labels);
    thrust::fill(dev_labels, dev_labels + graph.n, iterK);

    chkerr(cudaMemcpy(deviceGraph.degree, graph.degree.data(), graph.n * sizeof(ui), cudaMemcpyHostToDevice));

    //Tested


    size_t sharedMemoryIntialClique =  WARPS_EACH_BLK * sizeof(ui);
    listIntialCliques<<<1, 64,sharedMemoryIntialClique>>>(deviceDAG, levelData, labels, iterK, graph.n, graph.m, pSize, cpSize, maxBitMask, level, TOTAL_WARPS);
    CUDA_CHECK_ERROR("Generate Intial Partial Cliques");

    iterK--;
    level++;
    ui offsetPartitionSize = ((pSize / (k-1)) + 1);

    createLevelDataOffset(levelData, offsetPartitionSize, TOTAL_WARPS);
    
    flushParitions<<<BLK_NUMS, BLK_DIM>>>(*deviceDAG, *levelData, pSize, cpSize, maxBitMask, level, TOTAL_WARPS);
    CUDA_CHECK_ERROR("Flush Partition data structure");

    int totalTasks;
    chkerr(cudaMemcpy(&totalTasks, levelData->count + TOTAL_WARPS, sizeof(ui), cudaMemcpyDeviceToHost));

    while(iterK > 2) {
        thrust::fill(dev_labels, dev_labels + G.n, iterK);
        listMidCliques<<<BLK_NUMS, BLK_DIM>>>(*deviceDAG, *levelData, labels, k, iterK, G.n, G.m, pSize, cpSize, maxBitMask, totalTasks, level, TOTAL_WARPS);
        CUDA_CHECK_ERROR("Generate Mid Partial Cliques");

        createLevelDataOffset(*levelData, offsetPartitionSize, TOTAL_WARPS);

        flushParitions<<<BLK_NUMS, BLK_DIM>>>(*deviceDAG, *levelData, pSize, cpSize, k, maxBitMask, level, TOTAL_WARPS);
        CUDA_CHECK_ERROR("Flush Partition data structure");

        chkerr(cudaMemcpy(&totalTasks, levelData->count + TOTAL_WARPS, sizeof(ui), cudaMemcpyDeviceToHost));
        iterK--;
        level++;
    }

    ui t;
    //TODO: decide the total number cliques
    chkerr(cudaFree(labels));
    freeLevelPartitionData(*levelData);
    
    memoryAllocationTrie(*cliqueData, t, k);
    int totalCliques;

    chkerr(cudaMalloc((void**)&totalCliques, sizeof(ui)));
    chkerr(cudaMemset(totalCliques, 0, sizeof(ui)));

    if(iterK == 2) {
        writeFinalCliques<<<BLK_NUMS, BLK_DIM>>>(*deviceGraph, *deviceDAG, *levelData, *cliqueData, totalCliques, k, iterK, G.n, G.m, pSize, cpSize, maxBitMask, totalTasks, level, TOTAL_WARPS);
        CUDA_CHECK_ERROR("Generate Full Cliques");
    }

    freeLevelData(*levelData);

    sortTrieData<<<BLK_NUMS, BLK_DIM>>>(*deviceGraph, *cliqueData, totalCliques, k, TOTAL_THREAD);
    CUDA_CHECK_ERROR("Sort Trie Data Structure");

    freeDAG(*deviceDAG);

    // TODO: reorder Trie by motif degree

    //TODO:  CLIQUE CORE DECOMPOSE

    //TODO: LOCATE CORE

    //TODO: LISTING AGAIN

    //TODO: EDGE PRUNING

    //TODO: COMPONENT DECOMPOSE

    //TODO: DYNAMIC CORE EXACT

    freeGraph(*deviceGraph);
    delete deviceGraph;
    delete deviceDAG;
    delete levelData;
    delete cliqueData;
    return 0;
}