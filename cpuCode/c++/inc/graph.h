#pragma once

#include "common.h"

class Graph {
public:
  ui n;
  ui m;
  ui maxCliqueDegree;
  ui maxDegree;

  std::vector<std::vector<ui>> adjacencyList;
  std::vector<ui> degree;
  std::vector<ui> core;
  std::vector<long> cliqueDegree;
  std::vector<long> cliqueCore;
  std::string filePath;

public:
  Graph();
  Graph(std::string path);
};

class Motif {
public:
  ui size;
  ui type;
  ui count;
  std::string filePath;
  std::vector<std::vector<ui>> adjMatrix;

public:
  Motif();
  Motif(std::string path);
};