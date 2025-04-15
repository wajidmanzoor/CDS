#pragma once

#define miv(a, b) ((a) > (b) ? (b) : (a))
#define mav(a, b) ((a) < (b) ? (b) : (a))


#include <assert.h>
#include <string.h>

#include <cstdlib>
#include <fstream>
#include <list>
#include <queue>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <iostream>
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <set>
#include <map>
#include <utility>
#include <limits.h>
#include <ctime>
#include <sys/stat.h>
#include <mutex>

#define BLK_NUMS 1
#define BLK_DIM 64
#define TOTAL_THREAD (BLK_NUMS*BLK_DIM)
#define WARPSIZE 32
#define WARPS_EACH_BLK (BLK_DIM/32)
#define TOTAL_WARPS (BLK_NUMS*WARPS_EACH_BLK)

using namespace std;

typedef unsigned int ui;
typedef unsigned short ushort;
typedef unsigned char uchar;

const int INF = 1000000000;


typedef struct  {

    ui *offset;
    ui *neighbors;

    ui *degree;
    ui *cliqueDegree;
    
    int *cliqueCore;
    ui *cliqueCorePeelSequence;


    double *density;
    ui *motifCount;
}deviceGraphPointers;

typedef struct  {
    ui *adjacencyMatrix;
}deviceMotifPointers;

typedef struct  {
    ui *componentOffset;
    ui *components;
    ui *mapping;
}deviceComponentPointers;

typedef struct  {

    ui *maxDensity;
    ui *numVertex;
    ui *component;
    ui *status;
    
}deviceResultpointer;

typedef struct{
    ui *offset;
    ui *neighbors;
    ui *degree;
}deviceDAGpointer;

typedef struct {

    ui *trie;
    int *status;

}deviceCliquesPointer;

typedef struct {
    ui *partialCliquesPartition;
    ui  *partialCliques;         // Flattened array of partial cliques
    ui *candidatesPartition;       // Flattened array of candidate sets
    ui *candidates;       // Flattened array of candidate sets
    ui *offsetPartition;
    ui *offset; 
    ui *validNeighMaskPartition;   
    ui *validNeighMask;       
    ui *count; 
    ui *temp;           // Number of partial cliques at this level
    ui *max;
}cliqueLevelDataPointer; 

typedef struct{
    ui *mapping;
    ui *reverseMap;
    ui *offset;
    ui *neighbors;
    ui *cliqueDegree;
    double *density;
    ui *n;
    ui *m;
    ui *totalCliques;


}densestCorePointer;

typedef struct{
  ui *newOffset;
  ui *newNeighbors;
  ui *pruneStatus;

}devicePrunedNeighbors;

typedef struct{
    ui *toEdge;
    double *capacity;
    double *flow;

}deviceFlowNetworkPointers;

extern  deviceGraphPointers deviceGraph;
extern  deviceDAGpointer deviceDAG;
extern  cliqueLevelDataPointer levelData;
extern  deviceCliquesPointer cliqueData;
extern  densestCorePointer  densestCore;
extern  deviceComponentPointers conComp;
extern devicePrunedNeighbors prunedNeighbors;
extern deviceFlowNetworkPointers flowNetwork;
extern  deviceCliquesPointer finalCliqueData;




