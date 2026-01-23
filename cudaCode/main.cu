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

ui listAllCliquesBaseline(const Graph &graph, deviceDAGpointer &deviceDAG,
                          deviceCliquesPointer &cliqueData, ui k, ui pSize,
                          ui cSize) {

  cliqueLevelDataBaseline A, B;

  ui maxDegree = max_element(graph.degree.begin(), graph.degree.end());

  // TODO: CHECK
  ui maxBitMask = allocLevelDataBaseline(A, k, pSize, cSize, maxDegree);
  allocLevelDataBaseline(B, k, pSize, cSize, maxDegree);
  ui oneLabelSize = (graph.n + 31) / 32;

  size_t numWords = static_cast<size_t>(oneLabelSize) * TOTAL_WARPS;

  ui *labels;
  chkerr(cudaMalloc((void **)&(labels), numWords * sizeof(ui)));
  cudaMemset(labels, 0, numWords * sizeof(ui));

  chkerr(cudaMemcpy(deviceGraph.degree, graph.degree.data(),
                    graph.n * sizeof(ui), cudaMemcpyHostToDevice));

  // level 0
  listInitialCliquesBaseline<<<BLK_NUMS, BLK_DIM>>>(
      deviceDAG, A, label, k, graph.n, maxBitMask, TOTAL_WARPS);
  cudaDeviceSynchronize();

  cliqueLevelDataBaseline *read = &A, *write = &B;

  ui level = 1;

  ui totalTasks;

  while (level < k - 1) {

    listMidCliquesBaseline<<<BLK_NUMS, BLK_DIM>>>(D, *read, *write, label, k,
                                                  graph.n, level, maxBitMask,
                                                  level, TOTAL_WARPS);
    cudaDeviceSynchronize();

    std::swap(read, write);
    level++;
  }

  ui *globalCounter;
  cudaMalloc(&globalCounter, sizeof(ui));
  cudaMemset(globalCounter, 0, sizeof(ui));

  memoryAllocationTrie(cliqueData, MAX_TASKS, k);

  writeFinalCliquesBaseline<<<BLK_NUMS, BLK_DIM>>>(*read, deviceDAG, cliqueData,
                                                   globalCounter, k);
  cudaDeviceSynchronize();

  ui totalCliques;
  cudaMemcpy(&totalCliques, globalCounter, sizeof(ui), cudaMemcpyDeviceToHost);
  freeLevelDataBaseline(A);
  freeLevelDataBaseline(B);

  return totalCliques;
}

int main(int argc, const char *argv[]) {
  if (argc != 8) {
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
  ui glBufferSize =
      atoi(argv[5]); // Buffer to store the vertices that need to be removed
                     // in clique core decompose peeling algorithm
  ui partitionSize = atoi(argv[6]); // Virtual Partition to store the active
                                    // node of each flownetwork
  // ui t = atoi(argv[8]);             // Total Number of cliques.

  int earlyStop = atoi(argv[7]);

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
  ui coreTotalCliques, maxCore;
  double maxDensity;

  std::vector<ui> coreSize;
  std::chrono::high_resolution_clock::time_point coreDecomEnd;
  std::chrono::high_resolution_clock::time_point cliqueListEnd;
  std::chrono::duration<double, std::milli> cliqueList_ms;
  std::chrono::duration<double, std::milli> coreDecom_ms;

  totalCliques =
      listAllCliquesBaseline(graph, deviceDAG, cliqueData, k, pSize, cSize);

  cliqueListEnd = std::chrono::high_resolution_clock::now();
  cliqueList_ms = cliqueListEnd - dagEnd;

  freeGraph(deviceGraph);
  freeTrie(finalCliqueData);
  cout << "Total Cliques=" << totalCliques << endl;

  std::cout << "DAG Time =" << dag_ms.count() << " ms" << std::endl;
  std::cout << "Clique Listing Time =" << cliqueList_ms.count() << " ms"
            << std::endl;

  return 0;
}