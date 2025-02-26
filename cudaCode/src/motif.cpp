#include "../inc/motif.h"

Motif::Motif(){

}

Motif::motif(string filepath){
    path = filepath;
    string buffer;
    ifstream inputFile(path, uis::in);
    if(!inputFile.is_open()){
        cout << "Motif file Open Failed"<<endl;
        exit(1)
    }else{
        string line;
        getline(inputFile,line)
        istringstream iss(line);

        iss>> size >> type >> count;

        adjMatrix = new vector<ui>(size*size, 0);

        ui u,v;

        while(getline(inputFile,line)){
            istringstream iss(line);

            iss >> u >> v;

            adjMatrix[u*size][v] = 1;
            adjMatrix[v*size][u] = 1;
        }

    cout<< "Motif Size:"<<size<<", Motif Count: "<<count<<", Motif Type"<<type<<endl;
    }
}