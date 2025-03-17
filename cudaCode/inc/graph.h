#pragma once

#include "common.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <numeric>

class Graph {
  public:
    ui n;
    ui m;

    std::vector<ui> offset;
    std::vector<ui> neighbors;
    std::vector<ui> degree;
    std::vector<ui> core;
    std::vector<ui> corePeelSequence;
    std::string filePath;
  
  public:
  
    Graph();
    Graph(std::string path);
    void getListingOrder(std::vector<ui>& arr);
    void coreDecompose(std::vector<ui>& arr);
};

