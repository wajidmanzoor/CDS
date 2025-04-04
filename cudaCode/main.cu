#include "./inc/common.h"
#include "./inc/graph.h"

#include "./utils/cuda_utils.cuh"
#include "./inc/gpuMemoryAllocation.cuh"
#include "./inc/helpers.cuh"
#include <thrust/count.h>



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
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Generate Intial Partial Cliques");
    ui partialSize = TOTAL_WARPS * pSize;
    //ui candidateSize = TOTAL_WARPS * cpSize;
    ui offsetSize = ((pSize / (k - 1)) + 1) * TOTAL_WARPS;

  
    ui offsetPartitionSize = ((pSize / (k-1)) + 1);
    createLevelDataOffset(levelData, offsetPartitionSize, TOTAL_WARPS);

    flushParitions<<<BLK_NUMS, BLK_DIM>>>(deviceDAG, levelData, pSize,cpSize,k, maxBitMask, level, TOTAL_WARPS);
    cudaDeviceSynchronize();
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
        cudaDeviceSynchronize();
        CUDA_CHECK_ERROR("Generate Mid Partial Cliques");

        createLevelDataOffset(levelData, offsetPartitionSize, TOTAL_WARPS);

        chkerr(cudaMemset(levelData.offset,0,offsetSize*sizeof(ui)));
        chkerr(cudaMemset(levelData.validNeighMask,0,partialSize*maxBitMask*sizeof(ui)));

        flushParitions<<<BLK_NUMS, BLK_DIM>>>(deviceDAG, levelData, pSize,cpSize,k, maxBitMask, level, TOTAL_WARPS);
        cudaDeviceSynchronize();
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


    thrust::device_ptr<int> dev_ptr(cliqueData.status);

    thrust::fill(dev_ptr, dev_ptr + t, -2);
    chkerr(cudaMemset(cliqueData.trie, 0, t * k * sizeof(ui)));
    if(iterK == 2) {
        writeFinalCliques<<<BLK_NUMS, BLK_DIM,sharedMemoryFinal>>>(deviceGraph, deviceDAG, levelData, cliqueData, totalCliques, k, iterK, graph.n, graph.m, pSize, cpSize, maxBitMask, t,totalTasks, level, TOTAL_WARPS);
        cudaDeviceSynchronize();
        CUDA_CHECK_ERROR("Generate Full Cliques");
    }


    freeLevelData(levelData);
    freeLevelPartitionData(levelData);
    freeDAG(deviceDAG);

    ui tt;
    chkerr(cudaMemcpy(&tt, totalCliques, sizeof(ui), cudaMemcpyDeviceToHost));
    cout<<endl<<"total cliques "<<tt<<endl;

    // Sort Not needed
    //size_t sharedMemorySort =  2*k*WARPS_EACH_BLK * sizeof(ui);
    //sortTrieData<<<BLK_NUMS, BLK_DIM,sharedMemorySort>>>(deviceGraph, cliqueData, tt,t, k, TOTAL_THREAD);
    //cudaDeviceSynchronize();
    //CUDA_CHECK_ERROR("Sort Trie Data Structure");

    cudaFree(totalCliques);
    return tt;
}

void cliqueCoreDecompose(const Graph& graph,deviceGraphPointers& deviceGraph,deviceCliquesPointer& cliqueData, ui &maxCore, double &maxDensity, ui &coreSize, ui &coreTotalCliques, ui glBufferSize, ui k, ui t, ui tt){
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



    chkerr(cudaMemcpy(deviceGraph.cliqueCore, deviceGraph.cliqueDegree, graph.n * sizeof(ui), cudaMemcpyDeviceToDevice));


    thrust::device_vector<int> dev_vec(cliqueData.status, cliqueData.status + t);

    ui currentCliques = thrust::count(dev_vec.begin(), dev_vec.end(), -1);
    double currentDensity = static_cast<double>(currentCliques) / (graph.n - count);

    maxDensity = currentDensity;
    maxCore = 0;
    coreTotalCliques = currentCliques;
    coreSize = graph.n;


    while(count < graph.n){
        cudaMemset(bufTails, 0, sizeof(ui)*BLK_NUMS);

        // Select nodes whoes current degree is level, that means they should be removed as part of the level core
        selectNodes<<<BLK_NUMS, BLK_DIM>>>(deviceGraph, bufTails, glBuffers, glBufferSize, graph.n, level);
        cudaDeviceSynchronize();

        //Total number of verticies in buffer
        thrust::device_vector<ui> dev_vec1(bufTails, bufTails + BLK_NUMS);
        ui sum = thrust::reduce(dev_vec1.begin(), dev_vec1.end(), 0, thrust::plus<ui>());
        cout<<"sum "<<sum;
        cudaDeviceSynchronize();



        processNodesByWarp<<<BLK_NUMS, BLK_DIM>>>(deviceGraph, cliqueData , bufTails, glBuffers, globalCount, glBufferSize, graph.n, level, k, t, tt);
        cudaDeviceSynchronize();

        //debug
        chkerr(cudaMemcpy(h_coreCliques, deviceGraph.cliqueCore, graph.n * sizeof(int), cudaMemcpyDeviceToHost));
        cout<<endl<<"Core Cliques ";
        for(int i = 0; i < graph.n; i++){
            cout<<h_coreCliques[i]<<" ";
        }
        cout<<endl;
        //debug

        chkerr(cudaMemcpy(&count, globalCount, sizeof(unsigned int), cudaMemcpyDeviceToHost));

        level++;
        thrust::device_vector<int> dev_vec2(cliqueData.status, cliqueData.status + t);


       if((graph.n - count)!=0){
        currentCliques = thrust::count(dev_vec2.begin(), dev_vec2.end(), -1);
        currentDensity = static_cast<double>(currentCliques) / (graph.n - count);

        if(currentDensity>=maxDensity){
            maxDensity = currentDensity;
            maxCore = level;
            coreTotalCliques = currentCliques;
            coreSize = graph.n-count;

        }
        
    

       }

        cudaDeviceSynchronize();

    }
    cudaFree(globalCount);
    cudaFree(bufTails);
    cudaFree(glBuffers);

}


ui generateDensestCore(const Graph& graph,deviceGraphPointers& deviceGraph, densestCorePointer &densestCore, ui *reverseMap, ui coreSize, ui coreTotalCliques, ui lowerBoundDensity){
    memoryAllocationDensestCore(densestCore, coreSize, lowerBoundDensity , coreTotalCliques);
    ui *globalCount;

    chkerr(cudaMalloc((void**)&globalCount, sizeof(ui)));
    chkerr(cudaMemset(globalCount, 0, sizeof(ui)));


    generateDensestCore<<<BLK_NUMS, BLK_DIM>>>(deviceGraph,densestCore,globalCount,graph.n,lowerBoundDensity,TOTAL_WARPS);
    cudaDeviceSynchronize();

    thrust::inclusive_scan(thrust::device_ptr<ui>(densestCore.offset), thrust::device_ptr<ui>(densestCore.offset + coreSize + 1), thrust::device_ptr<ui>(densestCore.offset));

    ui edgeCountCore;
    chkerr(cudaMemcpy(&edgeCountCore, densestCore.offset+coreSize , sizeof(ui), cudaMemcpyDeviceToHost));
    chkerr(cudaMalloc((void**)&(densestCore.neighbors), edgeCountCore * sizeof(ui)));

    thrust::device_ptr<unsigned int> d_vertex_map_ptr(densestCore.mapping);

    thrust::device_ptr<unsigned int> d_reverse_map_ptr(reverseMap);

    thrust::device_vector<unsigned int> d_indices(coreSize);

    thrust::sequence(d_indices.begin(), d_indices.end());

    // Scatter indices into the reverse mapping array

    thrust::scatter(d_indices.begin(), d_indices.end(), d_vertex_map_ptr, d_reverse_map_ptr);

    size_t sharedMemoryGenNeighCore =  WARPS_EACH_BLK * sizeof(ui);
    generateNeighborDensestCore<<<BLK_NUMS, BLK_DIM,sharedMemoryGenNeighCore>>>(deviceGraph,densestCore,reverseMap,lowerBoundDensity,TOTAL_WARPS);
    cudaDeviceSynchronize();

    return edgeCountCore;

}

ui prune(densestCorePointer &densestCore, deviceCliquesPointer &cliqueData, ui *pruneStatus, ui *reverseMap, ui *newOffset, ui *newNeighbors, ui vertexCount, ui edgecount, ui k, ui t, ui t, ui lowerBoundDensity){
    
    //Prune
    chkerr(cudaMalloc((void**)&pruneStatus, edgecount * sizeof(ui)));

    thrust::device_ptr<ui> d_pruneStatus(pruneStatus);

    // Fill the array with 1 using Thrust
    thrust::fill(d_pruneStatus, d_pruneStatus + edgecount, 1);

    pruneEdges<<<BLK_NUMS, BLK_DIM>>>( densestCore,  cliqueData,reverseMap, pruneStatus, t, tt,  k, lowerBoundDensity);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Get prune status of each edge");


    // Get out degree after prune
    chkerr(cudaMalloc((void**)&newOffset, (vertexCount+1) * sizeof(ui)));
    chkerr(cudaMemset(newOffset, (vertexCount+1)  , sizeof(ui)));


    generateDegreeAfterPrune<<<BLK_NUMS, BLK_DIM>>>(densestCore , pruneStatus, newOffset, vertexCount , edgecount, TOTAL_WARPS);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Generate Degree after pruning");

    // cummulative sum offset
    thrust::inclusive_scan(thrust::device_ptr<ui>(newOffset), thrust::device_ptr<ui>(newOffset + vertexCount + 1), thrust::device_ptr<ui>(newOffset));

    ui newEdgeCount;
    chkerr(cudaMemcpy(&newEdgeCount, newOffset+ vertexCount , sizeof(ui), cudaMemcpyDeviceToHost));
    chkerr(cudaMalloc((void**)&(newNeighbors), newEdgeCount * sizeof(ui)));

    // Write neighbors of after output
    size_t sharedMemoryGenNeig =  WARPS_EACH_BLK * sizeof(ui);
    generateNeighborAfterPrune<<<BLK_NUMS, BLK_DIM,sharedMemoryGenNeig>>>(densestCore , pruneStatus, newOffset, newNeighbors, vertexCount, edgecount , TOTAL_WARPS);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Generate Neighbor after prune");

    return newEdgeCount;
    
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

    //Debug
    cout << "filepath: " << filepath << endl;
    cout << "motifPath: " << motifPath << endl;
    cout <<"k: " << k << endl;
    cout << "pSize: " << pSize << endl;
    cout << "cpSize: " << cpSize << endl;

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

    // CLIQUE CORE DECOMPOSE
    ui coreSize, coreTotalCliques,maxCore;
    double maxDensity;
    cout<<endl<<"BEFORE DECOMPOSE"<<endl;
    cliqueCoreDecompose(graph,deviceGraph,cliqueData,maxCore, maxDensity, coreSize, coreTotalCliques,glBufferSize, k,  t, tt);

    //LOCATE CORE
    ui *reverseMap;
    chkerr(cudaMalloc((void**)&reverseMap, graph.n * sizeof(ui)));
    cudaMemset(reverseMap, 0xFF, graph.n * sizeof(ui));
    ui lowerBoundDensity = static_cast<ui>(std::ceil(maxDensity));

    ui edgecount = generateDensestCore(graph,deviceGraph,  densestCore, reverseMap, coreSize, coreTotalCliques,lowerBoundDensity);

    //TODO: LISTING AGAIN not need added level as status of each clique to track the cores


    //TODO: EDGE PRUNING
    ui *pruneStatus;
    ui *newOffset;
    ui *newNeighbors;
    ui vertexCount;
    chkerr(cudaMemcpy(&vertexCount, densestCore.n, sizeof(ui), cudaMemcpyDeviceToHost));

    
    ui newEdgeCount = prune(densestCore, cliqueData, pruneStatus, reverseMap, newOffset, newNeighbors, vertexCount, edgecount, k, t, t, lowerBoundDensity);

    //TODO: COMPONENT DECOMPOSE

    //TODO: DYNAMIC CORE EXACT
    freTrie(cliqueData);
    freeGraph(deviceGraph);
    return 0;
}