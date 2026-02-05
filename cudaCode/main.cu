#include "./inc/common.h"
#include "./inc/graph.h"

#include "./inc/gpuMemoryAllocation.cuh"
#include "./inc/helpers.cuh"
#include "./utils/cuda_utils.cuh"
#include <thrust/count.h>

#include <cooperative_groups.h>
#include <cub/cub.cuh>
#include <iomanip>
#include <thrust/async/copy.h>
#include <thrust/binary_search.h>
#include <thrust/sequence.h>
#include <thrust/sort.h>
#include <thrust/unique.h>

#include <chrono>

bool DEBUG = false;
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <vector>

void generateDAG(const Graph &graph, deviceGraphPointers &deviceGraph,
                 deviceDAGpointer &deviceDAG, vector<ui> listingOrder) {

  // Stores the Directed Acyclic Graph
  memoryAllocationDAG(deviceDAG, graph.n, graph.m);

  ui *listOrder;
  chkerr(cudaMalloc((void **)&(listOrder), graph.n * sizeof(ui)));
  chkerr(cudaMemcpy(listOrder, listingOrder.data(), graph.n * sizeof(ui),
                    cudaMemcpyHostToDevice));

  // Get out degree in DAG
  generateDegreeDAG<<<BLK_NUMS, BLK_DIM>>>(deviceGraph, deviceDAG, listOrder,
                                           graph.n, graph.m, TOTAL_WARPS);
  cudaDeviceSynchronize();
  CUDA_CHECK_ERROR("Generate Degree of DAG");

  // copy out degree to offset
  chkerr(cudaMemset(deviceDAG.offset, 0, sizeof(ui)));
  chkerr(cudaMemcpy(deviceDAG.offset + 1, deviceDAG.degree,
                    (graph.n) * sizeof(ui), cudaMemcpyDeviceToDevice));

  // cummulative sum to get the offset of neighbors
  thrust::inclusive_scan(thrust::device_ptr<ui>(deviceDAG.offset),
                         thrust::device_ptr<ui>(deviceDAG.offset + graph.n + 1),
                         thrust::device_ptr<ui>(deviceDAG.offset));

  // Writes neighbors of DAG based on the offset
  size_t sharedMemoryGenDagNeig = WARPS_EACH_BLK * sizeof(ui);
  generateNeighborDAG<<<BLK_NUMS, BLK_DIM, sharedMemoryGenDagNeig>>>(
      deviceGraph, deviceDAG, listOrder, graph.n, graph.m, TOTAL_WARPS);
  cudaDeviceSynchronize();
  CUDA_CHECK_ERROR("Generate Neighbor of DAG");

  chkerr(cudaFree(listOrder));
}

ui listAllCliquesBaseline(const Graph &graph, deviceGraphPointers deviceGraph,
                          deviceDAGpointer &deviceDAG,
                          deviceCliquesPointer &cliqueData, ui k, ui pSize,
                          ui cSize) {

  cliqueLevelDataBaseline A, B;
  int iterK = k;

  chkerr(cudaMemcpy(deviceGraph.degree, graph.degree.data(),
                    graph.n * sizeof(ui), cudaMemcpyHostToDevice));

  thrust::device_ptr<ui> dev_degree(deviceDAG.degree);
  auto max_iter = thrust::max_element(dev_degree, dev_degree + graph.n);
  int maxDegree = *max_iter;

  // cout<<"MAX DEG: "<<maxDegree<<endl;

  // TODO: CHECK
  ui maxBitMask = allocLevelDataBaseline(A, k, pSize, cSize, maxDegree);

  cudaDeviceSynchronize();
  CUDA_CHECK_ERROR("Memory Allocation A");

  allocLevelDataBaseline(B, k, pSize, cSize, maxDegree);
  cudaDeviceSynchronize();
  CUDA_CHECK_ERROR("Memory Allocation B");
  ;
  ui oneLabelSize = (graph.n + 31) / 32;

  size_t numWords = static_cast<size_t>(oneLabelSize) * TOTAL_WARPS;

  ui *labels;
  chkerr(cudaMalloc((void **)&(labels), numWords * sizeof(ui)));
  cudaMemset(labels, 0, numWords * sizeof(ui));

  /*chkerr(cudaMemcpy(deviceGraph.degree, graph.degree.data(),
                    graph.n * sizeof(ui), cudaMemcpyHostToDevice));*/
  size_t sharedMemoryIntialClique = WARPS_EACH_BLK * sizeof(ui);

  // level 0

  ui *baseCounter;

  chkerr(cudaMalloc((void **)&baseCounter, sizeof(ui)));
  chkerr(cudaMemset(baseCounter, 0, sizeof(ui)));

  cudaDeviceSynchronize();

  listInitialCliquesBaseline<<<BLK_NUMS, BLK_DIM, sharedMemoryIntialClique>>>(
      deviceDAG, A, labels, k, graph.n, maxBitMask, TOTAL_WARPS, baseCounter);
  cudaDeviceSynchronize();
  CUDA_CHECK_ERROR("Generate Initial Partial Cliques");

  // debugPrintLevelDataHost_All_WithMaskCheck(A, deviceDAG, graph.n, k,
  // level=1, maxBitMask,
  // printBinaryMask=*/true);

  ui taskCountHost;
  cudaMemcpy(&taskCountHost, A.taskCount, sizeof(ui), cudaMemcpyDeviceToHost);

  ui lastOffset;
  cudaMemcpy(&lastOffset, &A.offset[taskCountHost], sizeof(ui),
             cudaMemcpyDeviceToHost);

  cudaMemcpy(&A.offset[taskCountHost + 1], &lastOffset, sizeof(ui),
             cudaMemcpyHostToDevice);

  cliqueLevelDataBaseline *read = &A, *write = &B;

  ui level = 1;
  iterK--;

  size_t sharedMemoryMid = WARPS_EACH_BLK * sizeof(ui);

  // cout<<"done with first one"<<endl;

  while (iterK > 2) {
    // cout<<"start with second one"<<endl;
    cudaMemset(write->taskCount, 0, sizeof(ui));
    cudaMemset(write->offset, 0,
               (pSize + 1) * sizeof(ui)); // at least offset[0]=0
    cudaMemset(write->lock, 0, sizeof(int));
    size_t mask_size = (size_t)cSize * maxBitMask * sizeof(ui);
    cudaMemset(write->validNeighMask, 0, mask_size);
    chkerr(cudaMemset(baseCounter, 0, sizeof(ui)));

    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Before  Mid Partial Cliques");

    cudaDeviceSynchronize();

    listMidCliquesBaseline<<<BLK_NUMS, BLK_DIM, sharedMemoryMid>>>(
        deviceDAG, *read, *write, labels, k, graph.n, maxBitMask, level,
        TOTAL_WARPS, baseCounter);

    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Generate Mid Partial Cliques");
    // cout<<"mid iter "<<iterK<<endl;

    // debugPrintLevelDataHost_All_WithMaskCheck(*write, deviceDAG, graph.n, k,
    //(level + 1), maxBitMask,
    // true);

    // Write offset[taskCount + 1] after final Mid level
    ui taskCountHost;
    cudaMemcpy(&taskCountHost, write->taskCount, sizeof(ui),
               cudaMemcpyDeviceToHost);

    ui lastOffset;
    cudaMemcpy(&lastOffset, &write->offset[taskCountHost], sizeof(ui),
               cudaMemcpyDeviceToHost);

    cudaMemcpy(&write->offset[taskCountHost + 1], &lastOffset, sizeof(ui),
               cudaMemcpyHostToDevice);
    std::swap(read, write);

    level++;
    iterK--;
  }

  ui *globalCounter;
  chkerr(cudaMalloc((void **)&globalCounter, sizeof(ui)));
  chkerr(cudaMemset(globalCounter, 0, sizeof(ui)));

  ui *totalCliques;
  chkerr(cudaMalloc((void **)&totalCliques, sizeof(ui)));
  chkerr(cudaMemset(totalCliques, 0, sizeof(ui)));
  cudaDeviceSynchronize();

  countCliques<<<BLK_NUMS, BLK_DIM>>>(deviceDAG, *read, totalCliques,
                                      maxBitMask, TOTAL_WARPS);
  cudaDeviceSynchronize();
  CUDA_CHECK_ERROR("Count Num Cliques");
  ui totalCliquesHost;
  cudaMemcpy(&totalCliquesHost, totalCliques, sizeof(ui),
             cudaMemcpyDeviceToHost);

  if (totalCliquesHost > 0) {
    memoryAllocationTrie(cliqueData, totalCliquesHost, k);

    // cout<<"TOTAL CLIQUES BEFOEE"<<totalCliquesHost<<endl;

    size_t sharedMemoryFinal = WARPS_EACH_BLK * sizeof(ui);

    writeFinalCliquesBaseline<<<BLK_NUMS, BLK_DIM, sharedMemoryFinal>>>(
        deviceGraph, *read, deviceDAG, cliqueData, globalCounter, k, maxBitMask,
        totalCliquesHost, TOTAL_WARPS);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Write Final Cliques");
  }

  freeLevelDataBaseline(A);
  freeLevelDataBaseline(B);

  return totalCliquesHost;
}

int main(int argc, const char *argv[]) {
  if (argc != 5) {
    cout << "Server wrong input parameters!" << endl;
    exit(1);
  }
  string filepath =
      argv[1]; //  Path to the graph file. The graph should be represented as
               //  an adjacency list with space separators
  // string motifPath = argv[2]; // Path to motif file. The motif should be
  // represented as edge list with space sperators. Not used yet
  ui k = atoi(argv[2]);      // The clique we are intrested in.
  ui pSize = atoi(argv[3]);  // Virtual Partition size for storing partial
                             // cliques in Listing Algorithm
  ui cpSize = atoi(argv[4]); // Virtual Partition size for storing Candidates
                             // of PC in Listing Algorithm
  // Read Graph as a adcajency List.
  Graph graph = Graph(filepath);

  // Stores the listing order based on core values:
  // vertices with higher core values are assigned lower (better) ranks.
  vector<ui> listingOrder;
  listingOrder.resize(graph.n);
  graph.getListingOrder(listingOrder);

  // Structure to store the graph on device
  memoryAllocationGraph(deviceGraph, graph);

  // Generates the DAG based on the listing order.
  // Only includes edges from a vertex with a lower listing order to one with
  // a higher listing order.
  auto start = std::chrono::high_resolution_clock::now();
  generateDAG(graph, deviceGraph, deviceDAG, listingOrder);
  auto dagEnd = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double, std::milli> dag_ms = dagEnd - start;

  ui totalCliques;
  // ui coreTotalCliques, maxCore;
  // double maxDensity;

  std::vector<ui> coreSize;
  std::chrono::high_resolution_clock::time_point cliqueListEnd;
  std::chrono::duration<double, std::milli> cliqueList_ms;
  totalCliques = listAllCliquesBaseline(graph, deviceGraph, deviceDAG,
                                        cliqueData, k, pSize, cpSize);

  cliqueListEnd = std::chrono::high_resolution_clock::now();
  cliqueList_ms = cliqueListEnd - dagEnd;

  /* ui cd[totalCliques*k];

   cudaMemcpy(cd,cliqueData.trie,totalCliques*k*sizeof(ui),cudaMemcpyDeviceToHost);

   cout<<"Clique data "<<endl;
   for(ui i=0;i<totalCliques;i++){
     for(ui j=0; j<k;j++){
       cout<<cd[j*totalCliques+i]<<" ";
     }
     cout<<endl;
   }*/

  freeGraph(deviceGraph);
  freeTrie(cliqueData);
  freeDAG(deviceDAG);
  // freeTrie(finalCliqueData);
  cout << "K: " << k << endl;
  cout << "Total Cliques=" << totalCliques << endl;

  std::cout << "DAG Time =" << dag_ms.count() << " ms" << std::endl;
  std::cout << "Clique Listing Time =" << cliqueList_ms.count() << " ms"
            << std::endl;

  return 0;
}