
__global__ void generateDegreeDAG(deviceGraphPointers &G, deviceDAGpointer &D, ui *listingOrder, ui n, ui m, ui totalWarps){


    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(ui i = warpId; i < n; i+=totalWarps){
        ui start = G.offset[i];
        ui end = G.offset[i+1];
        ui total = end - start;
        ui neigh;
        int count = 0;
        for(j = laneId; j < total; j +=warpSize){
            neigh = G.neighbors[start+j];
            if(listingOrder[i] < listingOrder[neigh]){
                count++;
            }

        }
        unsigned mask = __ballot_sync(0xFFFFFFFF, count > 0);
        int total_count = __popc(mask);
        
        if(laneId==0){
            D.degree[i] = total_count;
        }


    }
}

__global__ void generateNeighborDAG(deviceGraphPointers &G, deviceDAGpointer &D, ui *listingOrder, ui n, ui m, ui totalWarps){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(ui i = warpId; i < n; i+=totalWarps){
        ui start = G.offset[i];
        ui end = G.offset[i+1];
        ui total = end - start;
        ui neigh;
        int offset = G.offset[i];
        for(j = laneId; j < total; j +=warpSize){
            neigh = G.neighbors[start+j];
            if(listingOrder[i] < listingOrder[neigh]){
                int loc = atomicAdd(&offset,1);
                G.neighbors[loc] = neigh;
            }

        }
        
    }
}


