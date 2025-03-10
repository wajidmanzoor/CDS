#pragma once

#include "common.h"
#include "graph.h"
#include "motif.h"

// Function declarations for memory allocation
void memoryAllocationGraph(deviceGraphPointers &G, Graph &graph);
void memoryAllocationDAG(deviceDAGpointer &D, ui n, ui m);
void memoryAllocationMotif(deviceMotifPointers &M, Motif &motif);
void memoryAllocationComponent(deviceComponentPointers &C, ui n, ui m);
void memoryAllocationresult(deviceResultpointer &R, ui n);
void memeoryAllocationTrie(deviceCliquesPointer &C, ui t, ui k);
ui memoryAllocationlevelData(cliqueLevelDataPointer &L, ui k, ui pSize, ui cpSize, ui maxDegree, ui totalWarps);

// Function declarations for memory deallocation
void freeGraph(deviceGraphPointers &G);
void freeMotif(deviceMotifPointers &M);
void freeComponents(deviceComponentPointers &C);
void freeResults(deviceResultpointer &R);
void freTrie(deviceCliquesPointer &C);
void freeDAG(deviceDAGpointer &D);
void freeLevelPartitionData(cliqueLevelDataPointer &L);
void freeLevelData(cliqueLevelDataPointer &L);
