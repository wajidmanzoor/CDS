#include "../inc/motif.h"

Motif::Motif(){

}

Motif::Motif(string filepath) {
    path = filepath;
    string buffer;
    ifstream inputFile(path, ios::in);
    if(!inputFile.is_open()) {
        cout << "Motif file Open Failed" << endl;
        exit(1);
    } else {
        string line;
        getline(inputFile, line);
        istringstream iss(line);

        iss >> size >> type >> count;

        adjMatrix = new vector<vector<ui>>(size, vector<ui>(size, 0));

        ui u, v;

        while(getline(inputFile, line)) {
            istringstream iss(line);

            iss >> u >> v;

            (*adjMatrix)[u][v] = 1;
            (*adjMatrix)[v][u] = 1;
        }

        cout << "Motif Size:" << size << ", Motif Count: " << count << ", Motif Type" << type << endl;
    }
}