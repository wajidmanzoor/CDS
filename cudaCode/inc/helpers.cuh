#pragma once

#include "common.h"
#include "graph.h"

// Kernel function declarations
__global__ void generateDegreeDAG(deviceGraphPointers G, deviceDAGpointer D,
                                  ui *listingOrder, ui n, ui m, ui totalWarps);
__global__ void generateNeighborDAG(deviceGraphPointers G, deviceDAGpointer D,
                                    ui *listingOrder, ui n, ui m,
                                    ui totalWarps);
__global__ void listInitialCliquesBaseline(deviceDAGpointer D,
                                           cliqueLevelDataBaseline levelData,
                                           ui *label, ui k, ui n, ui maxBitMask,
                                           ui totalWarps);
__global__ void listMidCliquesBaseline(deviceDAGpointer D,
                                       cliqueLevelDataBaseline levelDataRead,
                                       cliqueLevelDataBaseline levelDataWrite,
                                       ui *label, ui k, ui n, ui maxBitMask,
                                       ui level, ui totalWarps);
__global__ void writeFinalCliquesBaseline(cliqueLevelDataBaseline levelData,
                                          deviceDAGpointer D,
                                          deviceCliquesPointer C,
                                          ui *globalCounter, ui k);
