#pragma once

#include "common.h"
#include "graph.h"

// Kernel function declarations
__global__ void generateDegreeDAG(deviceGraphPointers G, deviceDAGpointer D, ui *listingOrder, ui n, ui m, ui totalWarps);
__global__ void generateNeighborDAG(deviceGraphPointers G, deviceDAGpointer D, ui *listingOrder, ui n, ui m, ui totalWarps);
__global__ void listIntialCliques(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui *label, ui k, ui n, ui m, ui psize, ui cpSize, ui maxBitMask, ui level, ui totalWarps);
__global__ void flushParitions(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui psize, ui cpSize, ui k, ui maxBitMask, ui level, ui totalWarps);
__global__ void listMidCliques(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui *label, ui k, ui iterK, ui n, ui m, ui psize, ui cpSize, ui maxBitMask, ui totalTasks, ui level, ui totalWarps);
__global__ void writeFinalCliques(deviceGraphPointers G, deviceDAGpointer D, cliqueLevelDataPointer levelData, deviceCliquesPointer cliqueData, ui *globalCounter, ui k, ui iterK, ui n, ui m, ui psize, ui cpSize, ui maxBitMask, ui trieSize,ui totalTasks, ui level, ui totalWarps);
__global__ void sortTrieData(deviceGraphPointers G, deviceCliquesPointer cliqueData, ui totalCliques, ui t, ui k, ui totalThreads)
