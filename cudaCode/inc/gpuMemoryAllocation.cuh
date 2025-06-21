#pragma once

#include "common.h"
#include "graph.h"

// Function declarations for memory allocation
void memoryAllocationGraph(deviceGraphPointers &G, Graph &graph);
void memoryAllocationDAG(deviceDAGpointer &D, ui n, ui m);
// void memoryAllocationMotif(deviceMotifPointers &M, Motif &motif);
void memoryAllocationComponent(deviceComponentPointers &C, ui n, ui m);
// void memoryAllocationresult(deviceResultpointer &R, ui n);
void memoryAllocationTrie(deviceCliquesPointer &C, ui t, ui k);
ui memoryAllocationlevelData(cliqueLevelDataPointer &L, ui k, ui pSize,
                             ui cpSize, ui maxDegree, ui totalWarps);
void memoryAllocationDensestCore(densestCorePointer &C, ui n, ui density,
                                 ui totalCliques, ui graphsize);

void memoryAllocationPrunnedNeighbors(devicePrunedNeighbors &prunedNeighbors,
                                      ui n, ui m);

void memoryAllocationFlowNetwork(deviceFlowNetworkPointers &flowNetwork,
                                 ui vertexSize, ui neighborSize,
                                 ui totalComponents);

// Function declarations for memory deallocation
void freeGraph(deviceGraphPointers &G);
// void freeMotif(deviceMotifPointers &M);
void freeComponents(deviceComponentPointers &C);
// void freeResults(deviceResultpointer &R);
void freTrie(deviceCliquesPointer &C);
void freeDAG(deviceDAGpointer &D);
void freeLevelData(cliqueLevelDataPointer &L);
void freeDensestCore(densestCorePointer &C);
void freePruneneighbors(devicePrunedNeighbors &P);
void freeFlownetwork(deviceFlowNetworkPointers &F);
