#include "../utils/cuda_utils.cuh"
#include "../inc/gpuMemoryAllocation.cuh"


void memoryAllocationGraph(deviceGraphPointers &G, Graph &graph) {
    ui n = graph.n;
    ui m = graph.m;
    chkerr(cudaMalloc((void**)&(G.offset), (n + 1) * sizeof(ui)));
    chkerr(cudaMemcpy(G.offset, graph.offset.data(), (n + 1) * sizeof(ui), cudaMemcpyHostToDevice));

    chkerr(cudaMalloc((void**)&(G.neighbors), (2 * m) * sizeof(ui)));
    chkerr(cudaMemcpy(G.neighbors, graph.neighbors.data(), (2 * m) * sizeof(ui), cudaMemcpyHostToDevice));

    chkerr(cudaMalloc((void**)&(G.degree), n * sizeof(ui)));
    chkerr(cudaMemcpy(G.degree, graph.degree.data(), n * sizeof(ui), cudaMemcpyHostToDevice));

    chkerr(cudaMalloc((void**)&(G.cliqueDegree), n * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(G.cliqueCore), n * sizeof(int)));

    chkerr(cudaMalloc((void**)&(G.cliqueCorePeelSequence), n * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(G.density), n * sizeof(double)));
    chkerr(cudaMalloc((void**)&(G.motifCount), n * sizeof(ui)));

    cudaDeviceSynchronize();
}

void memoryAllocationDAG(deviceDAGpointer &D, ui n, ui m) {
    chkerr(cudaMalloc((void**)&(D.offset), (n + 1) * sizeof(ui)));
    chkerr(cudaMemset(D.offset, 0, (n + 1) * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(D.neighbors), m * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(D.degree), n * sizeof(ui)));
    cudaDeviceSynchronize();
}

/*void memoryAllocationMotif(deviceMotifPointers &M, Motif &motif) {
    ui n = motif.size;

    chkerr(cudaMalloc((void**)&(M.adjacencyMatrix), (n * n) * sizeof(ui)));

    // TODO: Flatten motif.adjMatrix
    std::vector<ui> flatMatrix(n * n);
    // Flatten the 2D adjacency matrix into a 1D array
    for (ui i = 0; i < n; ++i) {
        for (ui j = 0; j < n; ++j) {
            flatMatrix[i * n + j] = (*motif.adjMatrix)[i][j];
        }
    }

    chkerr(cudaMemcpy(M.adjacencyMatrix, flatMatrix.data(), (n * n) * sizeof(ui), cudaMemcpyHostToDevice));
    cudaDeviceSynchronize();
}*/

void memoryAllocationComponent(deviceComponentPointers &C,ui totalComponents, ui n, ui m) {
    chkerr(cudaMalloc((void**)&(C.componentOffset), (totalComponents + 1) * sizeof(ui)));
    chkerr(cudaMemset(C.componentOffset, 0, (totalComponents + 1) * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(C.components), n * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(C.mapping), n * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(C.reverseMapping), n * sizeof(ui)))
    cudaDeviceSynchronize();
}

void memoryAllocationresult(deviceResultpointer &R, ui n) {
    chkerr(cudaMalloc((void**)&(R.maxDensity), sizeof(double)));
    chkerr(cudaMalloc((void**)&(R.numVertex), sizeof(ui)));
    chkerr(cudaMalloc((void**)&(R.component), sizeof(ui)));
    chkerr(cudaMalloc((void**)&(R.status), n * sizeof(ui)));
    cudaDeviceSynchronize();
}

void memoryAllocationTrie(deviceCliquesPointer &C, ui t, ui k) {
    chkerr(cudaMalloc((void**)&(C.trie), (t * k) * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(C.status), t * sizeof(int)));
    cudaDeviceSynchronize();
}

ui memoryAllocationlevelData(cliqueLevelDataPointer &L, ui k, ui pSize, ui cpSize, ui maxDegree, ui totalWarps) {
    ui partialSize = totalWarps * pSize;
    ui candidateSize = totalWarps * cpSize;
    ui offsetSize = ((pSize / (k - 1)) + 1) * totalWarps;
    ui maxBitMask = (maxDegree + 31) / 32;
    ui maskSize = (cpSize * maxBitMask) * totalWarps;
    ui max_ = partialSize / (k - 1);

    chkerr(cudaMalloc((void**)&(L.partialCliquesPartition), partialSize * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(L.partialCliques), partialSize * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(L.candidatesPartition), candidateSize * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(L.candidates), candidateSize * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(L.validNeighMaskPartition), maskSize * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(L.validNeighMask), maskSize * sizeof(ui)));

    chkerr(cudaMemset(L.validNeighMask, 0, maskSize * sizeof(ui)));
    chkerr(cudaMemset(L.validNeighMaskPartition, 0, maskSize * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(L.offsetPartition), offsetSize * sizeof(ui)));
    chkerr(cudaMemset(L.offsetPartition, 0, offsetSize * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(L.offset), offsetSize * sizeof(ui)));
    chkerr(cudaMemset(L.offset, 0, offsetSize * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(L.count), (totalWarps + 1) * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(L.temp), (totalWarps + 1) * sizeof(ui)));
    chkerr(cudaMemset(L.temp, 0, (totalWarps + 1) * sizeof(ui)));
    chkerr(cudaMemset(L.count, 0, (totalWarps + 1) * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(L.max), sizeof(ui)));
    chkerr(cudaMemcpy(L.max, &max_, sizeof(ui), cudaMemcpyHostToDevice));
    cudaDeviceSynchronize();
    return maxBitMask;
}
void memoryAllocationDensestCore(densestCorePointer &C, ui n, ui density, ui totalCliques){


    chkerr(cudaMalloc((void**)&(C.mapping), n* sizeof(ui)));

    chkerr(cudaMalloc((void**)&(C.offset), (n+1) * sizeof(ui)));
    chkerr(cudaMemset(C.offset, 0, (n+1)* sizeof(ui)));

    //neighbors will be allocated once we now the size

    chkerr(cudaMalloc((void**)&(C.cliqueDegree), n * sizeof(ui)));
    //chkerr(cudaMalloc((void**)&(C.cliqueCore), n * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(C.density), n * sizeof(double)));
    chkerr(cudaMemcpy(C.density, &density, sizeof(double), cudaMemcpyHostToDevice));

    chkerr(cudaMalloc((void**)&(C.n), sizeof(ui)));
    chkerr(cudaMemcpy(C.n, &n, sizeof(ui), cudaMemcpyHostToDevice));
    chkerr(cudaMalloc((void**)&(C.m), sizeof(ui)));
    
    chkerr(cudaMalloc((void**)&(C.totalCliques), sizeof(ui)));
    chkerr(cudaMemcpy(C.totalCliques, &totalCliques, sizeof(ui), cudaMemcpyHostToDevice));  
    chkerr(cudaMalloc((void**)&(C.reverseMap), n * sizeof(ui)));
    cudaDeviceSynchronize();

}


// Memory deallocation functions
void freeGraph(deviceGraphPointers &G) {
    chkerr(cudaFree(G.offset));
    chkerr(cudaFree(G.neighbors));
    chkerr(cudaFree(G.degree));
    chkerr(cudaFree(G.cliqueDegree));
    chkerr(cudaFree(G.cliqueCore));
    chkerr(cudaFree(G.cliqueCorePeelSequence));
    chkerr(cudaFree(G.density));
    chkerr(cudaFree(G.motifCount));
}

void memoryAllocationPrunnedNeighbors(devicePrunedNeighbors &prunedNeighbors, ui n, ui m){
    chkerr(cudaMalloc((void**)&(prunedNeighbors.newOffset), (n+1) * sizeof(ui)));
    chkerr(cudaMemset(prunedNeighbors.newOffset, 0, (n+1)* sizeof(ui)));

    chkerr(cudaMalloc((void**)&(prunedNeighbors.pruneStatus), (2 * m) * sizeof(ui)));
    cudaDeviceSynchronize();
}



void memoryAllocationFlowNetwork(deviceFlowNetworkPointers &flowNetwork, ui vertexSize, ui neighborSize, ui totalComponents){

    chkerr(cudaMalloc((void**)&(flowNetwork.offset), (1 + totalComponents) * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(flowNetwork.neighborOffset1), (1 + totalComponents) * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(flowNetwork.neighborOffset2), (vertexSize) * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(flowNetwork.Edges), neighborSize * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(flowNetwork.capacityForward), neighborSize * sizeof(capacity)));
    chkerr(cudaMalloc((void**)&(flowNetwork.capacityBackward), neighborSize * sizeof(capacity)));

    chkerr(cudaMalloc((void**)&(flowNetwork.flowForward), neighborSize* sizeof(capacity)));
    chkerr(cudaMalloc((void**)&(flowNetwork.flowBackward), neighborSize* sizeof(capacity)));

    chkerr(cudaMalloc((void**)&(flowNetwork.height), (vertexSize-totalComponents)* sizeof(capacity)));
    chkerr(cudaMalloc((void**)&(flowNetwork.excess), (vertexSize-totalComponents)* sizeof(capacity)));
    cudaDeviceSynchronize();

}

/*void freeMotif(deviceMotifPointers &M) {
    chkerr(cudaFree(M.adjacencyMatrix));
}*/


void freeComponents(deviceComponentPointers &C) {
    chkerr(cudaFree(C.componentTotal));
    chkerr(cudaFree(C.componentOffset));
    chkerr(cudaFree(C.components));
    chkerr(cudaFree(C.mapping));

}

void freeResults(deviceResultpointer &R) {
    chkerr(cudaFree(R.maxDensity));
    chkerr(cudaFree(R.numVertex));
    chkerr(cudaFree(R.component));
    chkerr(cudaFree(R.status));
}

void freTrie(deviceCliquesPointer &C) {
    chkerr(cudaFree(C.trie));
    chkerr(cudaFree(C.status));
}

void freeDAG(deviceDAGpointer &D) {
    chkerr(cudaFree(D.offset));
    chkerr(cudaFree(D.neighbors));
    chkerr(cudaFree(D.degree));
}

void freeLevelPartitionData(cliqueLevelDataPointer &L) {
    chkerr(cudaFree(L.partialCliquesPartition));
    chkerr(cudaFree(L.candidatesPartition));
    chkerr(cudaFree(L.offsetPartition));
    chkerr(cudaFree(L.validNeighMaskPartition));
    chkerr(cudaFree(L.temp));
}

void freeLevelData(cliqueLevelDataPointer &L) {
    chkerr(cudaFree(L.partialCliques));
    chkerr(cudaFree(L.candidates));
    chkerr(cudaFree(L.offset));
    chkerr(cudaFree(L.validNeighMask));
    chkerr(cudaFree(L.count));
    chkerr(cudaFree(L.max));
}

void freeDensestCore(densestCorePointer &C){
    chkerr(cudaFree(C.mapping));
    chkerr(cudaFree(c.reverseMap));
    chkerr(cudaFree(C.offset));
    chkerr(cudaFree(C.neighbors));
    chkerr(cudaFree(C.density));
    chkerr(cudaFree(C.n));
    chkerr(cudaFree(C.m));
    chkerr(cudaFree(C.totalCliques));
    chkerr(cudaFree(C.cliqueDegree));
 
}

void freePrunnedNeighbors(devicePrunedNeighbors &prunedNeighbors){
    chkerr(cudaFree(prunedNeighbors.newOffset));
    chkerr(cudaFree(prunedNeighbors.newNeighbors));
    chkerr(cudaFree(prunedNeighbors.pruneStatus));
}

void freeFlowNetwork(deviceFlowNetworkPointers &flowNetwork){
    chkerr(cudaFree(flowNetwork.toEdge));
    chkerr(cudaFree(flowNetwork.capacity));
    chkerr(cudaFree(flowNetwork.flow));
    cudaDeviceSynchronize();

}

