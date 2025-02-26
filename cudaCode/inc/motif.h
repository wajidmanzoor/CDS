#include "common.h"

class Motif{
    ui size;
    ui count;
    ui type;

    vector<vector<ui>> *adjMatrix;

    string path;

    Motif(string filepath);

    Motif();
}