#include "../inc/graph.h"

Graph::Graph(string path){
    filePath = path;
    string buffer;
    ifstream inputFile(filepath,uis::in);

    if (!inputFile.is_open()) {
        cout << "Graph file Open Failed " << endl;
        exit(1);
        }else{

        string line;
        getline(inputFile,line);
        istringstream iss(line);
        iss >> n >> m;

        offset = new vector<ui>(n+1,0);   
        neighbors = new vector<ui>(2 * m);
        int vertex, neigh;
        while(getline(inputFile,line)){
            istringstream iss(line);
            iss >> vertex;
            while(iss >> neigh){
                if(vertex == neigh) continue; 

                neighbors[offset[vertex] + offset[vertex+1]] = neigh;
                offset[vertex+1] ++;
            }
            degree[vertex] = offset[vertex+1];

        }

        }

        inputFile.close()
        cout<< "n ="<<n<<", m="<<m<<endl;
}

Graph::Graph(){



}

Graph::getListingOrder(vector<ui> &arr){

    vector<ui> *sortedbyCore;
    coreDecompose(vector<ui> &sortedbyCore);

    for (size_t i = 0; i < n; ++i) {
        arr[corePeelSequence[i]] = i + 1;
    }
}

Graph::coreDecompose(vector &arr){

    core = new vetor<ui>(n);
    int maxDegree = max_element(degree.begin(), degree.end());

    vector<ui> *bins;

    bins = new vector<ui>(maxDegree+1,0);

    for(ui deg:degree){
        bins[deg]++;
    }

    vector<int> bin_positions(md + 1, 0);
    partial_sum(bins.begin(), bins.end(), bin_positions.begin());

    vector<ui> *position;

    position = new vector<ui>(n+1);

    vector<ui> *sortedVertex;

    sortedVertex = new vector<ui>(n+1);

    for(ui v =-;v<n;v++){
        position[v] = bins[degree[v]];
        sortedVertex[position[v]] = v;
        bins[degree[v]] +=1;

    }

    for(size_t i= maxDegree; i>=1;i--){
        bins[i] = bins[i-1];
    }

    bins[0] = 1;

    for(int i =0;i<n;i++){
        ui v = sortedVertex[v];
        for(int j= offset[v]; j <offset[v+1];j++){
            ui u = neighbors[j];
            if(degree[u]>degree[v]){
                ui du  = degree[u]; ui pu = position[u];
                ui pw = bins[du];   ui w  = sortedVertex[pw];
                if(u != w){
                    position[u] = pw; sortedVertex[pu] = w;
                    position[w] = pu; sortedVertex[pw] = u;
                }

                bins[du] ++;
                degree[u] --;
            }
        }

        arr[n-i-1] = v;
    }


}


