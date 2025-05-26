#include "./inc/common.h"
#include "./inc/graph.h"

#include "./utils/cuda_utils.cuh"
#include "./inc/gpuMemoryAllocation.cuh"
#include "./inc/helpers.cuh"
#include <thrust/count.h>

#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/unique.h>
#include <thrust/binary_search.h>
#include <thrust/async/copy.h>
#include <cooperative_groups.h>
#include <cub/cub.cuh>


bool DEBUG = true;

__global__ void printmap(densestCorePointer densestCore, ui coreSize){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if(idx<coreSize){
    printf("idx %d map %d \n",idx,densestCore.mapping[idx]);
  }
}

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


    if(DEBUG){
        ui *h_degree,*h_offset;
        h_degree = new ui[graph.n];
        h_offset = new ui[graph.n + 1];

        chkerr(cudaMemcpy(h_degree, deviceDAG.degree, graph.n * sizeof(ui), cudaMemcpyDeviceToHost));
        chkerr(cudaMemcpy(h_offset, deviceDAG.offset, (graph.n + 1) * sizeof(ui), cudaMemcpyDeviceToHost));

        cout<<endl<<endl<<"DAG DATA"<<endl;
        cout<<endl<<"DAG"<<endl<<"Degree ";
        for(int i = 0; i < graph.n; i++) {
            cout << h_degree[i] << " ";
        }
        cout<<endl<<"offset ";
        for(int i = 0; i < graph.n + 1; i++) {
            cout << h_offset[i] << " ";
        }
        cout<<endl;

    }
    

    // Write neighbors of DAG
    size_t sharedMemoryGenDagNeig =  WARPS_EACH_BLK * sizeof(ui);
    generateNeighborDAG<<<BLK_NUMS, BLK_DIM,sharedMemoryGenDagNeig>>>(deviceGraph, deviceDAG, listOrder, graph.n, graph.m, TOTAL_WARPS);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Generate Neighbor of DAG");

    if(DEBUG){
        ui *h_neighbors;
        h_neighbors = new ui[graph.m];
        chkerr(cudaMemcpy(h_neighbors, deviceDAG.neighbors, graph.m * sizeof(ui), cudaMemcpyDeviceToHost));
        cout<<"neigh ";
        for(int i = 0; i < graph.m; i++) {
            cout << h_neighbors[i] << " ";
        }
        cout<<endl;


    }
    


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

      
      chkerr(cudaMemcpy(&totalTasks, levelData.count + TOTAL_WARPS, sizeof(ui), cudaMemcpyDeviceToHost));

      iterK--;
      level++;
    }


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
        CUDA_CHECK_ERROR("Generate Full Cliques");
    }

    if(DEBUG){
        int *h_cliques,*status;
        h_cliques = new int[t*k];
        status = new int[t];
        chkerr(cudaMemcpy(h_cliques, cliqueData.trie, k * t * sizeof(ui), cudaMemcpyDeviceToHost));
        chkerr(cudaMemcpy(status, cliqueData.status, t * sizeof(ui), cudaMemcpyDeviceToHost));
        cout<<endl;

        cout<<endl<<"Cliques Data "<<endl;
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
        //ui v = 4;

        //chkerr(cudaMemcpy(deviceGraph.cliqueDegree + 7,&v ,  sizeof(ui), cudaMemcpyHostToDevice));
        chkerr(cudaMemcpy(h_cdegree, deviceGraph.cliqueDegree, graph.n* sizeof(ui), cudaMemcpyDeviceToHost));


        cout<<endl<<"Clique Degree"<<endl;

        for(int i = 0; i < graph.n; i++) {
            cout <<i<<" ";
        }
        cout<<endl;
        for(int i = 0; i < graph.n; i++) {
            cout << h_cdegree[i] << " ";
        }

        

        //size_t sharedMemorySort =  2*k*WARPS_EACH_BLK * sizeof(ui);
        //sortTrieData<<<BLK_NUMS, BLK_DIM,sharedMemorySort>>>(deviceGraph, cliqueData, tt,t, k, TOTAL_THREAD);
        //CUDA_CHECK_ERROR("Sort Trie Data Structure");


        cout<<endl;

    }
    ui tt;
    chkerr(cudaMemcpy(&tt, totalCliques, sizeof(ui), cudaMemcpyDeviceToHost));
    
    freeLevelData(levelData);
    freeLevelPartitionData(levelData);
    freeDAG(deviceDAG);

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
        cudaDeviceSynchronize();


        processNodesByWarp<<<BLK_NUMS, BLK_DIM>>>(deviceGraph, cliqueData , bufTails, glBuffers, globalCount, glBufferSize, graph.n, level, k, t, tt);
        cudaDeviceSynchronize();

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

ui generateDensestCore(const Graph& graph,deviceGraphPointers& deviceGraph, densestCorePointer &densestCore, ui coreSize, ui coreTotalCliques, ui lowerBoundDensity){
    memoryAllocationDensestCore(densestCore, coreSize, lowerBoundDensity , coreTotalCliques);

    ui *globalCount;

    chkerr(cudaMalloc((void**)&globalCount, sizeof(ui)));
    chkerr(cudaMemset(globalCount, 0, sizeof(ui)));

    generateDensestCore<<<BLK_NUMS, BLK_DIM>>>(deviceGraph,densestCore,globalCount,graph.n,lowerBoundDensity,TOTAL_WARPS);
    cudaDeviceSynchronize();

    ui *h_offset;
    h_offset = new ui[coreSize+1];

    chkerr(cudaMemcpy(h_offset, densestCore.offset, (coreSize+1) * sizeof(ui), cudaMemcpyDeviceToHost));

     cout<<endl<<"offset b ";
    for(int i=0;i<=coreSize;i++){
      cout<<h_offset[i]<<" ";
    }
    cout<<endl;



    thrust::inclusive_scan(thrust::device_ptr<ui>(densestCore.offset), thrust::device_ptr<ui>(densestCore.offset + coreSize + 1), thrust::device_ptr<ui>(densestCore.offset));

    //debug
    ui *h_mapping;
    h_mapping = new ui[coreSize];
    chkerr(cudaMemcpy(h_mapping, densestCore.mapping, coreSize * sizeof(ui), cudaMemcpyDeviceToHost));
    chkerr(cudaMemcpy(h_offset, densestCore.offset, (coreSize+1) * sizeof(ui), cudaMemcpyDeviceToHost));


    cout<<endl<<"Densest core data "<<endl;
    cout<<endl<<"mapping ";
    for(int i=0;i<coreSize;i++){
      cout<<h_mapping[i]<<" ";
    }
    cout<<endl<<"offset ";
    for(int i=0;i<=coreSize;i++){
      cout<<h_offset[i]<<" ";
    }
    cout<<endl;


    ui edgeCountCore;
    chkerr(cudaMemcpy(&edgeCountCore, densestCore.offset+coreSize , sizeof(ui), cudaMemcpyDeviceToHost));

    cout<<"edgeCountCore "<<edgeCountCore<<endl;

    chkerr(cudaMemcpy(densestCore.m,&edgeCountCore, sizeof(ui), cudaMemcpyHostToDevice));
    chkerr(cudaMalloc((void**)&(densestCore.neighbors), edgeCountCore * sizeof(ui)));

    thrust::device_ptr<unsigned int> d_vertex_map_ptr(densestCore.mapping);
    thrust::device_ptr<unsigned int> d_reverse_map_ptr(densestCore.reverseMap);


    thrust::device_vector<unsigned int> d_indices(coreSize);
    thrust::sequence(d_indices.begin(), d_indices.end());
    // Scatter indices into the reverse mapping array
    thrust::scatter(d_indices.begin(), d_indices.end(), d_vertex_map_ptr, d_reverse_map_ptr);


    cout<<"gaph size "<<graph.n<<" core size "<<coreSize<<endl;


    size_t sharedMemoryGenNeighCore =  WARPS_EACH_BLK * sizeof(ui);
    generateNeighborDensestCore<<<BLK_NUMS, BLK_DIM,sharedMemoryGenNeighCore>>>(deviceGraph,densestCore,lowerBoundDensity,TOTAL_WARPS);
    cudaDeviceSynchronize();

    //Debug
    ui *h_neighbors;
    h_neighbors = new ui[edgeCountCore];
    chkerr(cudaMemcpy(h_neighbors, densestCore.neighbors, edgeCountCore * sizeof(ui), cudaMemcpyDeviceToHost));

    cout<<endl<<"neighbors ";
    for(int i=0;i<edgeCountCore;i++){
      cout<<h_neighbors[i]<<" ";
    }
    cout<<endl;


    return edgeCountCore;


}
ui prune(densestCorePointer &densestCore, deviceCliquesPointer &cliqueData, devicePrunedNeighbors &prunedNeighbors,
    ui vertexCount, ui edgecount, ui k, ui t, ui tt, ui lowerBoundDensity) {

    // Allocate and initialize pruneStatus
    thrust::device_ptr<ui> d_pruneStatus(prunedNeighbors.pruneStatus);
    thrust::fill(d_pruneStatus, d_pruneStatus + edgecount, 1);

    // Kernel to determine pruning

    pruneEdges<<<BLK_NUMS, BLK_DIM>>>(densestCore, cliqueData,prunedNeighbors.pruneStatus, t, tt, k, lowerBoundDensity);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Get prune status of each edge");

    if(DEBUG){
        ui *h_pstatus = new ui[edgecount];
        chkerr(cudaMemcpy(h_pstatus, prunedNeighbors.pruneStatus, edgecount * sizeof(ui), cudaMemcpyDeviceToHost));

        cout<<endl<<"Neigh Status"<<endl;
        for (ui i = 0; i < edgecount; i++) {
            std::cout << h_pstatus[i] << " ";
        }
        std::cout << std::endl;

        delete[] h_pstatus;

    }
    


    // Allocate and initialize newOffset

    // Kernel to generate out-degrees
    generateDegreeAfterPrune<<<BLK_NUMS, BLK_DIM>>>(
        densestCore, prunedNeighbors.pruneStatus, prunedNeighbors.newOffset, vertexCount, edgecount, TOTAL_WARPS);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Generate Degree after pruning");

    // Inclusive scan to build offset array
    thrust::inclusive_scan(
        thrust::device_pointer_cast(prunedNeighbors.newOffset),
        thrust::device_pointer_cast(prunedNeighbors.newOffset + vertexCount + 1),
        thrust::device_pointer_cast(prunedNeighbors.newOffset));

    // Get total number of remaining edges
    ui newEdgeCount;
    chkerr(cudaMemcpy(&newEdgeCount, prunedNeighbors.newOffset + vertexCount, sizeof(ui), cudaMemcpyDeviceToHost));

    // Allocate memory for newNeighbors
    chkerr(cudaMalloc((void**)&(prunedNeighbors.newNeighbors), newEdgeCount * sizeof(ui)));

    // Kernel to generate neighbors list
    size_t sharedMemoryGenNeig = WARPS_EACH_BLK * sizeof(ui);
    generateNeighborAfterPrune<<<BLK_NUMS, BLK_DIM, sharedMemoryGenNeig>>>(
        densestCore, prunedNeighbors.pruneStatus, prunedNeighbors.newOffset, prunedNeighbors.newNeighbors, vertexCount, edgecount, TOTAL_WARPS);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Generate Neighbor after prune");


    if(DEBUG){
        ui *h_offset, *h_neigh;
        h_offset = new ui[vertexCount+1];
        h_neigh = new ui[newEdgeCount];
        chkerr(cudaMemcpy(h_offset, prunedNeighbors.newOffset, (vertexCount+1) * sizeof(ui), cudaMemcpyDeviceToHost));
        chkerr(cudaMemcpy(h_neigh, prunedNeighbors.newNeighbors, newEdgeCount * sizeof(ui), cudaMemcpyDeviceToHost));

        cout<<endl<<"Data After Pruning "<<endl;
        cout<<"in offset";
        for(ui i=0;i<vertexCount+1;i++){
        cout<<h_offset[i]<<" ";
        }cout<<endl;
        cout<<"in neigh";
        for(ui i=0;i<newEdgeCount;i++){
        cout<<h_neigh[i]<<" ";
        }cout<<endl;

    }
    

    return newEdgeCount;
}

int componentDecompose(deviceComponentPointers &conComp,devicePrunedNeighbors &prunedNeighbors, ui vertexCount, ui edgecount){

    ui *changed;
    chkerr(cudaMalloc((void**)&changed, sizeof(ui)));
    chkerr(cudaMemset(changed, 0 , sizeof(ui)));

    thrust::device_ptr<ui> components = thrust::device_pointer_cast(conComp.components);
    thrust::sequence(components, components + vertexCount);

    //int iter = 0;
    //can be used to put a limit on num iters
    ui hostChanged;
    do{
      chkerr(cudaMemset(changed, 0 , sizeof(ui)));
      componentDecomposek<<<BLK_NUMS, BLK_DIM>>>(conComp, prunedNeighbors, changed, vertexCount, edgecount, TOTAL_WARPS);
      cudaDeviceSynchronize();
      CUDA_CHECK_ERROR("Coponenet Decompose");

      chkerr(cudaMemcpy(&hostChanged,changed , sizeof(ui), cudaMemcpyDeviceToHost));

    }while(hostChanged>0);

    // unique component
    thrust::device_vector<int> uniqueComponents(vertexCount);
    auto new_end = thrust::unique_copy(components , components + vertexCount,
                                   uniqueComponents.begin());
    int totalComponents = new_end - uniqueComponents.begin();

    //Create component offset
    thrust::device_ptr<ui> componentOffsets(conComp.componentOffset);
    thrust::lower_bound(components , components + vertexCount,
                    uniqueComponents.begin(), uniqueComponents.begin() + totalComponents,
                    componentOffsets);
    componentOffsets[totalComponents] = vertexCount;

    uniqueComponents.resize(totalComponents);


    thrust::sort(uniqueComponents.begin(), uniqueComponents.end());


    thrust::lower_bound(
        uniqueComponents.begin(), uniqueComponents.end(),
        components, components + vertexCount,
        components // In-place remap
    );

    //Vertices in Densest Core, use mapping to get actual verticies
    thrust::device_ptr<ui> vertices(densestCore.mapping);

    thrust::sequence(thrust::device_pointer_cast(conComp.mapping), thrust::device_pointer_cast(conComp.mapping + vertexCount));

    thrust::sort_by_key(
      components,
      components + vertexCount,
      thrust::device_pointer_cast(conComp.mapping)
    );
    //TODO : new neighbor list, new neighbor offset

     return totalComponents;

}









int main(int argc, const char * argv[]) {
    if (argc != 8) {
        cout << "Server wrong input parameters!" << endl;
        exit(1);
    }

    string filepath = argv[1]; // Path to the graph file. The graph should be represented as an adjacency list with space separators
    string motifPath = argv[2]; //Path to motif file. The motif should be represented as edge list with space sperators.
    ui k = atoi(argv[3]);
    ui pSize = atoi(argv[4]);
    ui cpSize = atoi(argv[5]);
    ui glBufferSize = atoi(argv[6]);
    ui partitionSize = atoi(argv[7]);

    if(DEBUG){
        cout << "filepath: " << filepath << endl;
        cout << "motifPath: " << motifPath << endl;
        cout <<"k: " << k << endl;
        cout << "pSize: " << pSize << endl;
        cout << "cpSize: " << cpSize << endl;
    }
    //find a way to do this
    ui t = 10;
    
    Graph graph = Graph(filepath);

    //Print Graph
    if(DEBUG){
        cout<<"Graph Data "<<endl;
        cout<<"Graph"<<endl<<"Offset: ";
        for (int i = 0; i < (graph.n +1); i++) {
            cout<<graph.offset[i]<<" ";

        }
        cout<<endl<<"Neighbors: ";
        for (int i = 0; i < 2*graph.m; i++) {
            cout<<graph.neighbors[i]<<" ";
        }
        cout<<endl;
        cout<<"Degree: ";
        for (int i = 0; i < graph.n; i++) {
            cout<<graph.degree[i]<<" ";
        }
        cout<<endl;

    }
    
    vector<ui> listingOrder;
    listingOrder.resize(graph.n);
    graph.getListingOrder(listingOrder);

    if(DEBUG){
        cout<<endl<<endl<<"Listing Order: ";
        for (int i = 0; i < graph.n; i++) {
            cout << listingOrder[i] << " ";
        }
        cout<<endl<<"Core ";
        for(int i = 0; i < graph.n; i++) {
            cout << graph.core[i] << " ";
        }
        cout<<endl<<"Peel Seq ";
        for(int i = 0; i < graph.n; i++) {
            cout  << graph.corePeelSequence[i] << " ";
        }
        cout<<endl;

    }

    memoryAllocationGraph(deviceGraph, graph);

    generateDAG(graph,deviceGraph, deviceDAG,listingOrder);

    ui tt = listAllCliques(graph, deviceGraph, deviceDAG, levelData, k, pSize, cpSize,t);

    ui coreSize, coreTotalCliques,maxCore;
    double maxDensity;
    cliqueCoreDecompose(graph,deviceGraph,cliqueData,maxCore, maxDensity, coreSize, coreTotalCliques,glBufferSize, k,  t, tt);

    if(DEBUG){
        cout<<endl<<"Clique data after core decompose "<<endl;
        int *h_cliques,*status;
        h_cliques = new int[t*k];
        status = new int[t];
        chkerr(cudaMemcpy(h_cliques, cliqueData.trie, k * t * sizeof(ui), cudaMemcpyDeviceToHost));
        chkerr(cudaMemcpy(status, cliqueData.status, t * sizeof(ui), cudaMemcpyDeviceToHost));
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
        cout<<endl;
    }
    
    ui lowerBoundDensity = maxCore;
 
    ui edgecount = generateDensestCore(graph,deviceGraph,  densestCore, coreSize, coreTotalCliques,lowerBoundDensity);

    ui vertexCount;
    chkerr(cudaMemcpy(&vertexCount, densestCore.n, sizeof(ui), cudaMemcpyDeviceToHost));

    memoryAllocationPrunnedNeighbors(prunedNeighbors, vertexCount , edgecount);

    ui newEdgeCount = prune(densestCore, cliqueData, prunedNeighbors, vertexCount, edgecount, k, t, tt, lowerBoundDensity);

    memoryAllocationComponent(conComp, vertexCount , newEdgeCount);

    ui totalComponents = componentDecompose(conComp, prunedNeighbors, vertexCount, newEdgeCount);

    //Dynamic exact

    thrust::device_ptr<unsigned int> d_vertex_map_ptr(conComp.mapping);
    thrust::device_ptr<unsigned int> d_reverse_map_ptr(conComp.reverseMapping);

    thrust::device_vector<unsigned int> d_indices(vertexCount);
    thrust::sequence(d_indices.begin(), d_indices.end());
    thrust::scatter(d_indices.begin(), d_indices.end(), d_vertex_map_ptr, d_reverse_map_ptr);


    //Counter to store total cliques of each component
    ui *compCounter;
    chkerr(cudaMalloc((void**)&(compCounter), (totalComponents+1)* sizeof(ui)));
    chkerr(cudaMemset(compCounter, 0, (totalComponents+1) * sizeof(ui)));

    // Get total cliques of each connected components

    getConnectedComponentStatus<<<BLK_NUMS, BLK_DIM>>>(conComp,cliqueData, densestCore,compCounter,t, tt,k, maxCore,TOTAL_THREAD);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Calculate total cliques for each Component");

    thrust::inclusive_scan(
        thrust::device_pointer_cast(compCounter),
        thrust::device_pointer_cast(compCounter + totalComponents + 1),
        thrust::device_pointer_cast(compCounter));


    ui totaLCliques;
    chkerr(cudaMemcpy(&totaLCliques,compCounter+totalComponents,sizeof(ui),cudaMemcpyDeviceToHost));

    // Allocate memory for new clique data arranged my connected component
    memoryAllocationTrie(finalCliqueData, totaLCliques, k);


    ui *counter;
    chkerr(cudaMalloc((void**)&(counter), totalComponents* sizeof(ui)));
    chkerr(cudaMemset(counter, 0, totalComponents * sizeof(ui)));
    rearrangeCliqueData<<<BLK_NUMS, BLK_DIM>>>(conComp, cliqueData,  finalCliqueData, densestCore, compCounter, counter, t,  tt,  k, totaLCliques,TOTAL_THREAD);

    if(DEBUG){
        int *h_cliques,*status;
        h_cliques = new int[t*k];
        status = new int[t];
        cout<<endl<<"Clique data after rearange "<<endl;
        chkerr(cudaMemcpy(h_cliques, finalCliqueData.trie, k * totaLCliques * sizeof(ui), cudaMemcpyDeviceToHost));
        chkerr(cudaMemcpy(status, finalCliqueData.status, totaLCliques* sizeof(ui), cudaMemcpyDeviceToHost));

        for(int i =0;i<k;i++){
        cout<<endl<<"CL "<<i<<"  ";
        for(int j =0;j<totaLCliques;j++){
            cout<<h_cliques[i*totaLCliques+j]<<" ";
        }
        }
        cout<<endl<<"stat  ";
        for(int i = 0; i < totaLCliques; i++) {
            cout << status[i] << " ";

        }
        cout<<endl;

    }
    
    double* bounds;
    chkerr(cudaMalloc((void**)&bounds, (totalComponents*2)* sizeof(double)));
    double* upperBound = bounds;
    double* lowerBound = bounds + totalComponents;


    chkerr(cudaMalloc((void**)&(flowNetwork.offset), (1 + totalComponents) * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(flowNetwork.neighborOffset1), (1 + totalComponents) * sizeof(ui)));

    thrust::device_ptr<int> d_cliqueCore(deviceGraph.cliqueCore);

    // Find the maximum element using thrust::reduce
    int max_int = thrust::reduce(d_cliqueCore, d_cliqueCore + graph.n,
                               0, thrust::maximum<int>());

    // Convert to double and return
    double md =  static_cast<double>(max_int);

    getLbUbandSize<<<BLK_NUMS, BLK_DIM>>>( conComp, compCounter, lowerBound, upperBound, flowNetwork.offset, flowNetwork.neighborOffset1, totalComponents, k, md);

    thrust::inclusive_scan(
        thrust::device_pointer_cast(flowNetwork.offset),
        thrust::device_pointer_cast(flowNetwork.offset + totalComponents + 1),
        thrust::device_pointer_cast(flowNetwork.offset));
    thrust::inclusive_scan(
        thrust::device_pointer_cast(flowNetwork.neighborOffset1),
        thrust::device_pointer_cast(flowNetwork.neighborOffset1 + totalComponents + 1),
        thrust::device_pointer_cast(flowNetwork.neighborOffset1));

    ui flownetworkVertexSize, flownetworkNeighborSize;
    chkerr(cudaMemcpy(&flownetworkVertexSize, flowNetwork.offset+totalComponents, sizeof(ui), cudaMemcpyDeviceToHost));
    chkerr(cudaMemcpy(&flownetworkNeighborSize, flowNetwork.neighborOffset1+totalComponents, sizeof(ui), cudaMemcpyDeviceToHost));

    memoryAllocationFlowNetwork(flowNetwork, flownetworkVertexSize, flownetworkNeighborSize, totalComponents);

    if(DEBUG){
        double *lb, *ub;
        lb = new double[totalComponents];
        ub = new double[totalComponents];

        ui *h_ccoffset, *h_neighborSize;

        h_ccoffset = new ui[totalComponents+1];
        h_neighborSize = new ui[totalComponents+1];

        chkerr(cudaMemcpy(lb,lowerBound,totalComponents*sizeof(double),cudaMemcpyDeviceToHost));
        chkerr(cudaMemcpy(ub,upperBound,totalComponents*sizeof(double),cudaMemcpyDeviceToHost));
        chkerr(cudaMemcpy(h_ccoffset,flowNetwork.offset,(totalComponents+1)*sizeof(ui),cudaMemcpyDeviceToHost));
        chkerr(cudaMemcpy(h_neighborSize,flowNetwork.neighborOffset1,(totalComponents+1)*sizeof(ui),cudaMemcpyDeviceToHost));

        cout<<endl<<"FLOW NETWORK"<<endl;
        cout<<"Lower Bound ";
        for(ui i=0;i<totalComponents;i++){
        cout<<lb[i]<<" ";
        }
        cout<<endl;

        cout<<"Upper Bound ";
        for(ui i=0;i<totalComponents;i++){
        cout<<ub[i]<<" ";
        }
        cout<<endl;
        cout<<"Vertex Offset ";
        for(ui i=0;i<=totalComponents;i++){
        cout<<h_ccoffset[i]<<" ";
        }
        cout<<endl;
        cout<<"Neigh Size ";
        for(ui i=0;i<=totalComponents;i++){
        cout<<h_neighborSize[i]<<" ";
        }
        cout<<endl;
    }
    thrust::device_ptr<double> dev_lowerBound(lowerBound);
     // Find the maximum element in the lowerBound array
    thrust::device_ptr<double> max_iter = thrust::max_element(dev_lowerBound, dev_lowerBound + totalComponents);

    // Copy the result back to host if needed
    double max_lowerBound = *max_iter;

    //double just  = 1.0;

    createFlowNetworkOffset<<<BLK_NUMS, BLK_DIM>>>( deviceGraph, flowNetwork, conComp, densestCore, finalCliqueData,   compCounter, upperBound,TOTAL_WARPS, totalComponents, k, max_lowerBound, totaLCliques);


    thrust::inclusive_scan(
        thrust::device_pointer_cast(flowNetwork.neighborOffset2),
        thrust::device_pointer_cast(flowNetwork.neighborOffset2 + flownetworkVertexSize),
        thrust::device_pointer_cast(flowNetwork.neighborOffset2));

    if(DEBUG){
        ui *h_offset2;
        h_offset2 = new ui[flownetworkVertexSize];

        chkerr(cudaMemcpy(h_offset2,flowNetwork.neighborOffset2,(flownetworkVertexSize)*sizeof(ui),cudaMemcpyDeviceToHost));

        cout<<"Neighbor Offset ";
        for(ui i=0;i<(flownetworkVertexSize);i++){
        cout<<h_offset2[i]<<" ";
        }
        cout<<endl;
    }
    createFlowNetwork<<<BLK_NUMS, BLK_DIM>>>( flowNetwork,  conComp,  densestCore,  finalCliqueData, compCounter,upperBound , TOTAL_WARPS, totalComponents, k, max_lowerBound , totaLCliques);
    cudaDeviceSynchronize();

    if(DEBUG){

        ui *neigh_s;
        neigh_s = new ui[flownetworkNeighborSize];

        chkerr(cudaMemcpy(neigh_s,flowNetwork.Edges,(flownetworkNeighborSize)*sizeof(ui),cudaMemcpyDeviceToHost));

        cout<<"neighbor ";
        for(ui i=0;i<(flownetworkNeighborSize);i++){
        cout<<neigh_s[i]<<" ";
        }
        cout<<endl;
        double *f_capacity;

        f_capacity =  new double[flownetworkNeighborSize];


        chkerr(cudaMemcpy(f_capacity,flowNetwork.capacity,(flownetworkNeighborSize)*sizeof(double),cudaMemcpyDeviceToHost));

        cout<<"Capacity ";
        for(ui i=0; i<flownetworkNeighborSize; i++){
        cout<<f_capacity[i]<<" ";
        }
        cout<<endl;

        ui *height;
        double *excess;
        height = new ui[flownetworkVertexSize];
        excess = new double[flownetworkVertexSize];

        chkerr(cudaMemcpy(height,flowNetwork.height, flownetworkVertexSize*sizeof(ui),cudaMemcpyDeviceToHost));
        chkerr(cudaMemcpy(excess,flowNetwork.excess, flownetworkVertexSize*sizeof(double),cudaMemcpyDeviceToHost));

        for(ui i =0;i <flownetworkVertexSize; i++){
        cout<<height[i]<<" ";
        }cout<<endl;

        for(ui i =0; i<flownetworkVertexSize; i++){
        cout<<excess[i]<<" ";
        }cout<<endl;

    }

    ui *activeNodes;

    cudaMalloc((void**)&activeNodes, partitionSize*TOTAL_WARPS * sizeof(ui));
    chkerr(cudaMemset(flowNetwork.flow,0,(flownetworkNeighborSize)* sizeof(double)));

    size_t sharedMemorySize=  WARPS_EACH_BLK * sizeof(ui)+ WARPS_EACH_BLK *sizeof(ui);
    ui *componenetsLeft;
    ui left = 1;

    chkerr(cudaMalloc((void**)&componenetsLeft,sizeof(ui)));
    chkerr(cudaMemset(componenetsLeft,0,sizeof(ui)));


    double *densities;
    cudaMalloc((void**)&densities,totalComponents*sizeof(double));
     chkerr(cudaMemset(densities,0,totalComponents*sizeof(double)));

     //float *densityMax;

    //ui wajid =0;
    while(left!=0){
      chkerr(cudaMemset(componenetsLeft,0,sizeof(ui)));
      pushRelabel<<<BLK_NUMS, BLK_DIM,sharedMemorySize>>>( flowNetwork,  conComp,  densestCore,  finalCliqueData, compCounter,upperBound,lowerBound, activeNodes,componenetsLeft,TOTAL_WARPS, totalComponents, k,totaLCliques, partitionSize);
      cudaDeviceSynchronize();
      cudaMemcpy(&left,componenetsLeft,sizeof(ui),cudaMemcpyDeviceToHost);
      
      getResult<<<BLK_NUMS, BLK_DIM,sharedMemorySize>>>( flowNetwork,  conComp,  finalCliqueData, compCounter, upperBound, lowerBound, densities, TOTAL_WARPS, totalComponents, k, t);
      cudaDeviceSynchronize();
      chkerr(cudaMemset(flowNetwork.flow,0,(flownetworkNeighborSize)* sizeof(double)));


    }

    cudaDeviceSynchronize();

    freeGraph(deviceGraph);


    return 0;
}