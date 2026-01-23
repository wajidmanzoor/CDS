#pragma once

#include "common.h"
#include "graph.h"

// Function declarations for memory allocation
void memoryAllocationGraph(deviceGraphPointers &G, Graph &graph);
void memoryAllocationDAG(deviceDAGpointer &D, ui n, ui m);

void memoryAllocationTrie(deviceCliquesPointer &C, ui t, ui k);

ui allocLevelDataBaseline(cliqueLevelDataBaseline &L, ui k, ui pSize, ui cSize,
                          ui maxDegree);

// Function declarations for memory deallocation
void freeGraph(deviceGraphPointers &G);
void freeDAG(deviceDAGpointer &D);
void freeTrie(deviceCliquesPointer &C);
void freeLevelDataBaseline(cliqueLevelDataBaseline &L);
