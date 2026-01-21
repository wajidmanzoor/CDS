#pragma once

#define miv(a, b) ((a) > (b) ? (b) : (a))
#define mav(a, b) ((a) < (b) ? (b) : (a))

#include <assert.h>
#include <cfloat>
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
#include <bitset>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits.h>
#include <map>
#include <mutex>
#include <set>
#include <sys/stat.h>
#include <utility>

#define debug 0
using namespace std;
using Clock = std::chrono::high_resolution_clock;

typedef unsigned int ui;
typedef unsigned short ushort;
typedef unsigned char uchar;

typedef unsigned char byte;
typedef unsigned long long ull;

const int inf = 1e9;
const double DINF = 1e9;

struct DensestCoreData {
  vector<vector<ui>> graph;
  vector<ui> reverseMap;
  int lowerBound;
  int delVertexIndex;
  int delCliqueCount;
  double density;
  int maxCliqueCore;
};

struct ConnectedComponentData {
  vector<vector<ui>> graph;
  int size;
  vector<ui> reverseMap;
  unordered_map<string, vector<int>> cliqueData;
  vector<long> cliqueDegree;
  long totalCliques;
  double density;
};

struct finalResult {
  vector<ui> verticies;
  double density;
  ui size;
};

struct Edge {
  int to;
  int rev;   // index of reverse edge
  float cap; // remaining capacity
  float max; // max capacity
};

struct FlowNetwork {
  vector<vector<Edge>> G;
  vector<int> sinkEdgeIdx; // cache vertex → sink edge
  int source, sink;

  void init(int n) {
    G.clear();
    G.resize(n);
    sinkEdgeIdx.assign(n, -1);
  }

  inline void addEdge(int u, int v, float cap) {
    Edge a{v, (int)G[v].size(), cap, cap};
    Edge b{u, (int)G[u].size(), 0.0f, 0.0f};
    G[u].push_back(a);
    G[v].push_back(b);
  }
};
