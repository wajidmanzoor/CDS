#include "../inc/common.h"
#include ".../inc/graph.h"

void memoryAllocationGraph(deviceGraphPointers &G, Graph &graph){

    ui n = graph.n;
    ui m = graph.m;
    
    chkerr(cudaMalloc((void**)&(G.offset), (n+1) * sizeof(ui)));
    chkerr(cudaMemcpy(G.offset, graph.offset->data(), (n+1) * sizeof(ui), cudaMemcpyHostToDevice));

    chkerr(cudaMalloc((void**)&(G.neighbors), (2*m) * sizeof(ui)));
    chkerr(cudaMemcpy(G.neighbors, graph.neighbors->data(), (2*m) * sizeof(ui), cudaMemcpyHostToDevice));

    chkerr(cudaMalloc((void**)&(G.degree), n * sizeof(ui)));
    chkerr(cudaMemcpy(G.degree, graph.degree->data(), n * sizeof(ui), cudaMemcpyHostToDevice));

    chkerr(cudaMalloc((void**)&(G.cliqueDegree), n * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(G.cliqueCorePeelSequence), n * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(G.density), n * sizeof(double)));

    chkerr(cudaMalloc((void**)&(G.motifCount), n * sizeof(ui)));

    cudaDeviceSynchronize();

}

void memoryAllocationDAG(deviceDAGpointer &D, ui n, ui m){
    chkerr(cudaMalloc((void**)&(D.offset), (n+1) * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(D.neighbors), m * sizeof(ui)));
    chkerr(cudaMalloc((void**)&(D.degree), n * sizeof(ui)));
    cudaDeviceSynchronize();
    
}

void memoryAllocationMotif(deviceMotifPointers &M, Motif &motif, ui n){

    ui n = motif.size;

    chkerr(cudaMalloc((void**)&(M.adjacencyMatrix), (n*n) * sizeof(ui)));
    chkerr(cudaMemcpy(M.adjacencyMatrix,motif.adjMatrix->data() , (n*n) * sizeof(ui), cudaMemcpyHostToDevice));
    cudaDeviceSynchronize();

}

void memoryAllocationComponent(deviceComponentPointers &C, ui n, ui m,){

    chkerr(cudaMalloc((void**)&(C.componentOffset), (n+1) * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(C.offset), (n+1) * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(C.neighbors), (2*m) * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(C.cliqueDegree), n * sizeof(ui)));


    chkerr(cudaMalloc((void**)&(C.cliqueCore), n * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(C.density), n * sizeof(double)));

    chkerr(cudaMalloc((void**)&(C.motifCount), n * sizeof(ui)));

    chkerr(cudaMalloc((void**)&(C.cliqueCorePeelSequence), n * sizeof(ui)));
    cudaDeviceSynchronize();

}

void memoryAllocationresult(deviceResultpointer &R){
    chkerr(cudaMalloc((void**)&(R.maxDensity), sizeof(double)));
    chkerr(cudaMalloc((void**)&(R.numVertex), sizeof(ui)));
    chkerr(cudaMalloc((void**)&(R.component), sizeof(ui)));
    chkerr(cudaMalloc((void**)&(R.status), (n)*sizeof(ui)));
    cudaDeviceSynchronize();

}


void memeoryAllocationTrie(deviceCliquesPointer &C, ui n, ui k){
    chkerr(cudaMalloc((void**)&(C.trie), (n*k)*sizeof(ui)));
    chkerr(cudaMalloc((void**)&(C.status),n*sizeof(ui)));
    cudaDeviceSynchronize();

}



void freeGraph(deviceGraphPointers &G){

    chkerr(cudaFree(G.offset));
    chkerr(cudaFree(G.neighbors));
    chkerr(cudaFree(G.degree));
    chkerr(cudaFree(G.cliqueDegree));
    chkerr(cudaFree(G.cliqueCore));
    chkerr(cudaFree(G.cliqueCorePeelSequence));
    chkerr(cudaFree(G.density));
    chkerr(cudaFree(G.motifCount));


}

void freeMotif(deviceMotifPointers &M){

    chkerr(cudaFree(M.adjacencyMatrix));

}

void freeComponents(deviceComponentPointers &C){

    chkerr(cudaFree(C.componentOffset));
    chkerr(cudaFree(C.offset));
    chkerr(cudaFree(C.neighbors));
    chkerr(cudaFree(C.cliqueDegree));
    chkerr(cudaFree(C.cliqueCore));
    chkerr(cudaFree(C.density));
    chkerr(cudaFree(C.motifCount));
    chkerr(cudaFree(C.cliqueCorePeelSequence));


}

void freeResults(deviceResultpointer &R){

    chkerr(cudaFree(R.maxDensity));
    chkerr(cudaFree(R.numVertex));
    chkerr(cudaFree(R.component));
    chkerr(cudaFree(R.status));
}

void freTrie(deviceCliquesPointer &C){
    chkerr(cudaFree(C.trie));
    chkerr(cudaFree(C.status));
}

void freeDAG(deviceDAGpointer &D){
    chkerr(cudaFree(D.offset));
    chkerr(cudaFree(D.neighbors));
    chkerr(cudaFree(D.degree));
}