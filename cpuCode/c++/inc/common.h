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
#include <bitset>
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
#define debug 1
using namespace std;

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