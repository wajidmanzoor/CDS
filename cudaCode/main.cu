#include "./inc/graph.h"
#include "./inc/motif.h"

#include "./src/gpuMemoryAllocation.cc"
#include "./src/helpers.cc"


#define CUDA_CHECK_ERROR(kernelName) { \
    cudaError_t err = cudaGetLastError(); \
    if (err != cudaSuccess) { \
        printf("CUDA Error in kernel %s, file %s at line %d: %s\n", kernelName, __FILE__, __LINE__, cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
    } \
}

bool fileExists(const string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

void writeOrAppend(const string& filename, const string& data) {
    ofstream file;
    
    // Check if the file exists
    if (fileExists(filename)) {
        // Open the file in append mode if it exists
        file.open(filename, ios::app);
    } else {
        // Open the file in write mode if it doesn't exist
        file.open(filename);
    }
    
    if (file.is_open()) {
        file << data << endl;
        file.close();
    } else {
        cerr << "Unable to open the file." << endl;
    }
}

int main(int argc, const char * argv[]) {
    if (argc != 3) {
        cout << "Server wrong input parameters!" << endl;
        exit(1);
    }

    string filepath = argv[1]; // Path to the graph file. The graph should be represented as an adjacency list with space separators
    string motifPath = argv[2]; //Path to motif file. The motif should be represented as edge list with space sperators.

    Graph G = Graph(filepath);

    Motif M = Motif(motifPath);

    vector<ui> *listingOrder;
    
    G.getListingOrder(&listingOrder);
    memoryAllocationGraph(deviceGraph,G);
    memoryAllocationDAG(deviceDAG,G.n;G.m);

    // Get out degree in DAG
    generateDegreeDAG<<<BLK_NUMS, BLK_DIM>>>(deviceGraph,deviceDAG,listingOrder,G.n, G.m,TOTAL_WARPS);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Generate Degree of DAG");

    //copy out degree to offset
    chkerr(cudaMemset(deviceDAG.neighbor, 0, sizeof(ui)));
    chkerr(cudaMemcpy(deviceDAG.offset + 1,deviceDAG.degree , (G.n) * sizeof(ui), cudaMemcpyDevicetoDevice));

    // cummulative sum offset
    thrust::inclusive_scan(thrust::device_ptr < ui > (deviceDAG.offset), thrust::device_ptr < ui > (deviceDAG.offset + G.n + 1), thrust::device_ptr < ui > (deviceDAG.offset ));

    // Write neighbors of DAG
    generateNeighborDAG<<<BLK_NUMS, BLK_DIM>>>(deviceGraph,deviceDAG,listingOrder,G.n, G.m,TOTAL_WARPS);
    cudaDeviceSynchronize();
    CUDA_CHECK_ERROR("Generate Neighbors of DAG");


    freeGraph(deviceGraph);
    freeDAG(deviceDAG);
    delete G;
    delete M;
    return 0;
}
