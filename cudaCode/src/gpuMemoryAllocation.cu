
#include "../inc/gpuMemoryAllocation.cuh"
#include "../utils/cuda_utils.cuh"

void memoryAllocationGraph(deviceGraphPointers &G, Graph &graph) {
  ui n = graph.n;
  ui m = graph.m;
  chkerr(cudaMalloc((void **)&(G.offset), (n + 1) * sizeof(ui)));
  chkerr(cudaMemcpy(G.offset, graph.offset.data(), (n + 1) * sizeof(ui),
                    cudaMemcpyHostToDevice));

  chkerr(cudaMalloc((void **)&(G.neighbors), (2 * m) * sizeof(ui)));
  chkerr(cudaMemcpy(G.neighbors, graph.neighbors.data(), (2 * m) * sizeof(ui),
                    cudaMemcpyHostToDevice));

  chkerr(cudaMalloc((void **)&(G.degree), n * sizeof(ui)));
  chkerr(cudaMemcpy(G.degree, graph.degree.data(), n * sizeof(ui),
                    cudaMemcpyHostToDevice));

  chkerr(cudaMalloc((void **)&(G.cliqueDegree), n * sizeof(ui)));
  chkerr(cudaMalloc((void **)&(G.cliqueCore), n * sizeof(int)));
  chkerr(cudaMalloc((void **)&(G.cliqueCorePeelSequence), n * sizeof(ui)));
  chkerr(cudaMalloc((void **)&(G.density), n * sizeof(double)));
  chkerr(cudaMalloc((void **)&(G.motifCount), n * sizeof(ui)));

  cudaDeviceSynchronize();
}

void memoryAllocationDAG(deviceDAGpointer &D, ui n, ui m) {
  chkerr(cudaMalloc((void **)&(D.offset), (n + 1) * sizeof(ui)));
  chkerr(cudaMemset(D.offset, 0, (n + 1) * sizeof(ui)));

  chkerr(cudaMalloc((void **)&(D.neighbors), m * sizeof(ui)));
  chkerr(cudaMalloc((void **)&(D.degree), n * sizeof(ui)));
  cudaDeviceSynchronize();
}

void memoryAllocationTrie(deviceCliquesPointer &C, ui t, ui k) {
  chkerr(cudaMalloc((void **)&(C.trie), (t * k) * sizeof(ui)));
  chkerr(cudaMalloc((void **)&(C.status), t * sizeof(int)));
  cudaDeviceSynchronize();
}

ui allocLevelDataBaseline(cliqueLevelDataBaseline &L, ui k, ui pSize, ui cSize,
                          ui maxDegree) {
  ui maxBitMask = (maxDegree + 31) / 32;

  cudaMalloc(&L.partialCliques, pSize * sizeof(ui));
  cudaMalloc(&L.candidates, cSize * sizeof(ui));
  cudaMalloc(&L.offset, cSize * sizeof(ui));
  cudaMalloc(&L.validNeighMask, cSize * maxBitMask * sizeof(ui));
  cudaMalloc(&L.taskCount, sizeof(ui));
  cudaMalloc(&L.lock, sizeof(ui));

  cudaMemset(L.taskCount, 0, sizeof(ui));
  cudaMemset(L.lock, 0, sizeof(ui));

  cudaMemset(L.offset, 0, cSize * sizeof(ui));

  return maxBitMask;
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

void freeTrie(deviceCliquesPointer &C) {
  chkerr(cudaFree(C.trie));
  chkerr(cudaFree(C.status));
}

void freeDAG(deviceDAGpointer &D) {
  chkerr(cudaFree(D.offset));
  chkerr(cudaFree(D.neighbors));
  chkerr(cudaFree(D.degree));
}

void freeLevelDataBaseline(cliqueLevelDataBaseline &L) {
  chkerr(cudaFree(L.partialCliques));
  chkerr(cudaFree(L.candidates));
  chkerr(cudaFree(L.offset));
  chkerr(cudaFree(L.validNeighMask));
  chkerr(cudaFree(L.taskCount));
}
