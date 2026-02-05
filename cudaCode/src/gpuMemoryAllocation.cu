
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

  // Calculate and log expected memory usage
  size_t total_bytes = 0;

  // Clear any previous CUDA errors
  cudaGetLastError();

  // 1. Allocate partialCliques
  cudaError_t err =
      cudaMalloc((void **)&(L.partialCliques), pSize * (k - 1) * sizeof(ui));
  if (err != cudaSuccess) {
    cout << "Failed to allocate partialCliques: " << cudaGetErrorString(err)
         << " (size: " << pSize * (k - 1) * sizeof(ui) / (1024 * 1024) << " MB)"
         << endl;
    return 0;
  }
  total_bytes += pSize * (k - 1) * sizeof(ui);

  // 2. Allocate candidates
  err = cudaMalloc((void **)&(L.candidates), cSize * sizeof(ui));
  if (err != cudaSuccess) {
    cout << "Failed to allocate candidates: " << cudaGetErrorString(err)
         << " (size: " << cSize * sizeof(ui) / (1024 * 1024) << " MB)" << endl;
    cudaFree(L.partialCliques);
    return 0;
  }
  total_bytes += cSize * sizeof(ui);

  // 3. Allocate offset - FIXED: should be (pSize + 1) based on your kernel
  // usage
  err = cudaMalloc((void **)&(L.offset), (pSize + 1) * sizeof(ui));
  if (err != cudaSuccess) {
    cout << "Failed to allocate offset: " << cudaGetErrorString(err)
         << " (size: " << (pSize + 1) * sizeof(ui) / (1024 * 1024) << " MB)"
         << endl;
    cudaFree(L.partialCliques);
    cudaFree(L.candidates);
    return 0;
  }
  total_bytes += (pSize + 1) * sizeof(ui);
  size_t mask_size = (size_t)cSize * maxBitMask * sizeof(ui);

  err = cudaMalloc((void **)&(L.validNeighMask), mask_size);
  if (err != cudaSuccess) {
    cout << "Failed to allocate validNeighMask: " << cudaGetErrorString(err)
         << " (size: " << mask_size / (1024 * 1024) << " MB)" << endl;
    cudaFree(L.partialCliques);
    cudaFree(L.candidates);
    cudaFree(L.offset);
    return 0;
  }
  total_bytes += mask_size;

  // 5. Allocate taskCount
  err = cudaMalloc((void **)&(L.taskCount), sizeof(ui));
  if (err != cudaSuccess) {
    cout << "Failed to allocate taskCount: " << cudaGetErrorString(err) << endl;
    cudaFree(L.partialCliques);
    cudaFree(L.candidates);
    cudaFree(L.offset);
    cudaFree(L.validNeighMask);
    return 0;
  }
  total_bytes += sizeof(ui);

  // 6. Allocate lock
  err = cudaMalloc((void **)&(L.lock), sizeof(int));
  if (err != cudaSuccess) {
    cout << "Failed to allocate lock: " << cudaGetErrorString(err) << endl;
    cudaFree(L.partialCliques);
    cudaFree(L.candidates);
    cudaFree(L.offset);
    cudaFree(L.validNeighMask);
    cudaFree(L.taskCount);
    return 0;
  }
  total_bytes += sizeof(int);

  // Initialize memory
  cudaMemset(L.taskCount, 0, sizeof(ui));
  cudaMemset(L.lock, 0, sizeof(int));
  cudaMemset(L.offset, 0, (pSize + 1) * sizeof(ui));
  cudaMemset(L.validNeighMask, 0, mask_size);

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
