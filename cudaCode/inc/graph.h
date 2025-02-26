
#pragma once

#include "common.h"

class Graph{

    ui n;
    ui m;

    vector<ui> *offset;
    vector<ui> *neighbors;
    vector<ui> *degree;
    vector<ui> *core;
    vector<ui> *corePeelSequence;

    string filePath;

    Graph(string path);

    Graph();

    void getListingOrder(vector<ui> &arr);

    void coreDecompose(vector<ui> &arr);

    

    
}
