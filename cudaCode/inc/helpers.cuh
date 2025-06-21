#pragma once

#include "common.h"
#include "graph.h"

// Kernel function declarations
__global__ void generateDegreeDAG(deviceGraphPointers G, deviceDAGpointer D, ui *listingOrder, ui n, ui m, ui totalWarps);
__global__ void generateNeighborDAG(deviceGraphPointers G, deviceDAGpointer D, ui *listingOrder, ui n, ui m, ui totalWarps);
__global__ void listIntialCliques(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui *label, ui k, ui n, ui m, ui psize, ui cpSize, ui maxBitMask, ui level, ui totalWarps);
__global__ void flushParitions(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui pSize, ui cpSize, ui k, ui maxBitMask, ui level, ui totalWarps);
__global__ void listMidCliques(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui *label, ui k,ui iterK, ui n, ui m,ui pSize, ui cpSize, ui maxBitMask,ui totalTasks, ui level, ui totalWarps);
__global__ void writeFinalCliques(deviceGraphPointers G, deviceDAGpointer D, cliqueLevelDataPointer levelData, deviceCliquesPointer cliqueData, ui *globalCounter,ui k,ui iterK, ui n, ui m,ui pSize, ui cpSize, ui maxBitMask,ui trieSize,ui totalTasks, ui level, ui totalWarps);
__global__ void countCliques(deviceDAGpointer D, cliqueLevelDataPointer levelData,  ui *globalCounter, ui maxBitMask ,ui totalTasks, ui totalWarps);

__global__ void sortTrieData(deviceGraphPointers G, deviceCliquesPointer cliqueData, ui totalCliques,ui t, ui k, ui totalThreads);

// Kernel function to Clique Core Decomposition
__global__ void selectNodes(deviceGraphPointers G, ui *bufTails, ui *glBuffers, ui glBufferSize, ui n, ui level);
__global__ void processNodesByWarp(deviceGraphPointers G,deviceCliquesPointer cliqueData, ui *bufTails,ui *glBuffers, ui *globalCount, ui glBufferSize, ui n, ui level, ui k, ui totalCliques);
__global__ void processNodesByBlock(deviceGraphPointers G,deviceCliquesPointer cliqueData, ui *bufTails,ui *glBuffers, ui *globalCount, ui glBufferSize, ui n, ui level, ui k, ui t, ui tt);


//kernel function to locate the densest core
__global__ void generateDensestCore(deviceGraphPointers G, densestCorePointer densestCore,ui *globalCount, ui n, ui core, ui totalWarps);
__global__ void generateNeighborDensestCore(deviceGraphPointers G, densestCorePointer densestCore, ui core, ui totalWarps);

__global__ void pruneEdges(densestCorePointer densestCore, deviceCliquesPointer cliqueData, ui *pruneStatus,ui totalCliques, ui k, ui level );
__global__ void generateDegreeAfterPrune(densestCorePointer densestCore ,ui *pruneStatus, ui *newOffset, ui n, ui m, ui totalWarps);
__global__ void generateNeighborAfterPrune(densestCorePointer densestCore ,ui *pruneStatus, ui *newOffset, ui *newNeighbors,ui n, ui m, ui totalWarps);

__global__ void componentDecomposek(deviceComponentPointers conComp, devicePrunedNeighbors prunedNeighbors, ui *changed, ui n, ui m, ui totalWarps);


//Dynamic exact
__global__ void getConnectedComponentStatus(deviceComponentPointers conComp,deviceCliquesPointer cliqueData, densestCorePointer densestCore, ui *compCounter, ui totalCliques, ui k,ui maxCore, ui totalThreads);
__global__ void rearrangeCliqueData(deviceComponentPointers conComp,deviceCliquesPointer cliqueData, deviceCliquesPointer finalCliqueData,densestCorePointer densestCore, ui *compCounter,ui *counter,ui totaLCliques,ui k, ui newTotaLCliques ,ui totalThreads);
__global__ void getLbUbandSize(deviceComponentPointers conComp, ui *compCounter, double *lowerBound, double *upperBound, ui *ccOffset,  ui *neighborSize, ui totalComponenets, ui k, double maxDensity);
__global__ void createFlowNetworkOffset(deviceGraphPointers G, deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp, deviceCliquesPointer finalCliqueData, ui *compCounter,double *upperBound , ui totalWarps, ui totalComponents, ui k, double lb, ui t);
__global__ void createFlowNetwork(deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp, deviceCliquesPointer finalCliqueData, ui *compCounter,double *upperBound , ui totalWarps, ui totalComponents, ui k, double lb, ui t);
__global__ void pushRelabel(deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp,  deviceCliquesPointer finalCliqueData, ui * compCounter, double * upperBound, double * lowerBound, ui * activeNodes, ui * componenetsLeft, ui * checkResult, ui totalWarps, int totalComponents, ui k, ui t, ui partitionSize);
__global__ void getResult(deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp, deviceCliquesPointer finalCliqueData, ui *compCounter,  double * upperBound, double * lowerBound, double *densities, ui *checkResult, ui totalWarps, int totalComponents, ui k,ui t);
 //__global__ void pushRelabel(deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp, densestCorePointer densestCore, deviceCliquesPointer finalCliqueData, ui *compCounter,double *upperBound, ui totalWarps, int totalComponents, ui k, double lb);

