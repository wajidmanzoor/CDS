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

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iostream>
#include <limits.h>
#include <map>
#include <mutex>
#include <set>
#include <sys/stat.h>
#include <utility>

#define NUM_BLKS_PER_SM 1
#define BLK_NUMS 216
#define BLK_DIM 1024
#define TOTAL_THREAD (BLK_NUMS * BLK_DIM)
#define WARPSIZE 32
#define WARPS_EACH_BLK (BLK_DIM / 32)
#define TOTAL_WARPS (BLK_NUMS * WARPS_EACH_BLK)

using namespace std;

typedef unsigned int ui;
typedef unsigned short ushort;
typedef unsigned char uchar;

const int INF = 1000000000;
const double DINF = 1e9;

typedef struct {

  ui *offset;
  ui *neighbors;

  ui *degree;
  ui *cliqueDegree;

  int *cliqueCore;
  ui *cliqueCorePeelSequence;

  double *density;
  ui *motifCount;
} deviceGraphPointers;

typedef struct {
  ui *offset;
  ui *neighbors;
  ui *degree;
} deviceDAGpointer;

typedef struct {

  ui *trie;
  int *status;

} deviceCliquesPointer;

typedef struct {
  ui *partialCliques; // [maxTasks * (k-1)]
  ui *candidates;     // [maxCandidates]
  ui *offset;         // [maxTasks + 1]
  ui *validNeighMask; // [maxCandidates * maxBitMask]

  ui *taskCount; // scalar on device
  ui *lock;
} cliqueLevelDataBaseline;

extern deviceGraphPointers deviceGraph;
extern deviceDAGpointer deviceDAG;
extern deviceCliquesPointer cliqueData;
