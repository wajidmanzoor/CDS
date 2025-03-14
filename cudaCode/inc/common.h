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


#define BLK_NUMS 432
#define BLK_DIM 512
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
    
    ui *cliqueCore;
    ui *cliqueCorePeelSequence;


    double *density;
    ui *motifCount;
}deviceGraphPointers;

typedef struct  {
    ui *adjacencyMatrix;
}deviceMotifPointers;

typedef struct  {

    ui *componentOffset;
    ui *offset;
    ui *neighbors;
    
    ui *cliqueDegree;
    ui *cliqueCore;
    
    double *density;
    ui *motifCount;
    ui *cliqueCorePeelSequence;
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
    ui *status;

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

extern  deviceGraphPointers deviceGraph;
extern  deviceDAGpointer deviceDAG;
extern  cliqueLevelDataPointer levelData;
extern  deviceCliquesPointer cliqueData;



