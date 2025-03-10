#pragma once
#include "common.h"

class Motif{

    public:
        ui size;
        ui count;
        ui type;
        string path;

    public:
        vector<vector<ui>> *adjMatrix;


    Motif(string filepath);
    Motif();
}