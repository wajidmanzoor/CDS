

__global__ void generateDegreeDAG(deviceGraphPointers G, deviceDAGpointer D, ui listingOrder, ui n, ui m, ui totalWarps){


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

__global__ void generateNeighborDAG(deviceGraphPointers G, deviceDAGpointer D, ui listingOrder, ui n, ui m, ui totalWarps){
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

__global__ void listIntialCliques(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui &label,ui k, ui n, ui m,ui psize, ui cpSize, ui maxBitMask, ui level, ui totalWarps){

    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui *counter = (ui * )(sharedMemory + sizeOffset);
   
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;
    int cliquePartition  = warpId*pSize;
    int offsetPartition = warpId*(pSize/(k-1)+1);
    int candidatePartition = warpId*cpSize;
    int maskPartition = WarpId*cpSize*maxBitMask; 

    
    for(int i =warpId; i < n ; i+= totalWarps ){
        ui vertex = i;
        ui neighOffset = D.offset[vertex];
        if(laneId==0){
            counter[threadId.x/warpSize] = 0;
        }

        __syncwarp();
        int candidateOffset =  candidatePartition + levelData.offsetPartition[offsetPartition + levelData.count[warpId+1]];
        for(int j = laneId; j< D.degree[vertex]; j+= warpSize ){
            ui neigh = D.neighbors[neighOffset + j];
            if(label[neigh] == k){
                label[neigh] = k-1;
                ui loc = atomicAdd( &counter[threadId.x/warpSize], 1);
                levelData.candidatesPartition[writeOffset + loc] = neigh;

            }

        }
        __syncwarp();
        if(laneId == 0 && counter[threadId.x/warpSize] > 0){
            levelData.partialCliquesPartition[cliquePartition + levelData.count[warpId+1] * (k-1) + level ] = vertex;
            levelData.count[warpId+1] +=1;
            levelData.offsetPartition[offsetPartition + levelData.count[warpId+1]] =
                levelData.offsetPartition[offsetPartition + levelData.count[warpId+1] - 1] +counter[threadId.x/warpSize];
        }
    }

    __syncwarp();

    int totalTasks = levelData.count[warpId+1];

    for(int iter = 0; iter < totalTasks; iter++){
        int start = candidatePartition + levelData.offsetPartition[offsetPartition + iter];
        int end =  candidatePartition + levelData.offsetPartition[offsetPartition + iter + 1];
        int total = end-start;
        for(int i = laneId; i < total; i+=warpSize){
            int candidate = levelData.candidatesPartition[start+i]; 
            int neighOffset = D.offset[candidate];
            int degree = D.degree[candidate];

            int numBitmasks = (degree + 31) / 32;

            for (int bitmaskIndex = 0; bitmaskIndex < numBitmasks; bitmaskIndex++) {
                unsigned int bitmask = 0; // Initialize bitmask to 0

                // Iterate over the current chunk of 32 neighbors
                int startNeighbor = bitmaskIndex * 32;
                int endNeighbor = min(startNeighbor + 32, degree);
                for (int j = startNeighbor; j < endNeighbor; j++) {
                    if (label[D.neighbors[neighOffset + j]] == k - 1) {
                        bitmask |= (1 << (j - startNeighbor)); // Set the bit for valid neighbors
                    }
                }

                D.validNeighMaskPartition[ maskPartition + (levelData.offsetPartition[offsetPartition + iter] + i)*maxBitMask + bitmaskIndex] = bitmask;
            
            }


        }

    }

}

__global__ void flushParitions(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui psize, ui cpSize, ui k, ui maxBitMask, ui level, ui totalWarps){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;
    int cliquePartition  = warpId*pSize;
    int offsetPartition = warpId*(pSize/(k-1)+1);
    int candidatePartition = warpId*cpSize;
    int maskPartition = WarpId*cpSize*maxBitMask; 

    int totalTasks = levelData.count[warpId+1] - levelData.count[warpId];

    for(int iter = 0; iter < totalTasks; iter++){
        int start = candidatePartition + levelData.offsetPartition[offsetPartition + iter];
        int end = candidatePartition + levelData.offsetPartition[offsetPartition + iter+ 1];
        int total end-start;
        int writeOffset = levelData.temp[warpId] + levelData.offsetPartition[offsetPartition + iter];
        for(int i = laneId; i < total; i+=warpSize){
            levelData.candidates[writeOffset+ i] = levelData.candidatesPartition[start + i];
            if(i< level){
                levelData.partialCliques[levelData.count[warpId]+ iter*(k-1) + i] = leveldata.partialCliquesPartition[cliquePartition + iter * (k-1) + i]
            }
            int totalMasks = (D.degree[candidate]+31)/32;
            for(int j =0; j < totalMasks; j++){
                levelData.validNeighMask[(writeOffset+i)*maxBitMask + j ] = 
                levelData.validNeighMaskPartition[maskPartition + (levelData.offsetPartition[offsetPartition + iter] + i)*maxBitMask + j]
            }

        }

        __syncWarp();

        if(laneId==0){
            levelData.offset[levelData.count[warpId] + iter] = levelData.temp[warpId]+levelData.offsetPartition[offsetPartition + iter+ 1];

        }


    }

}


__global__ void listMidCliques(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui &label,ui k,ui iterK, ui n, ui m,ui psize, ui cpSize, ui maxBitMask,ui totalTasks, ui level, ui totalWarps){

    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui *counter = (ui * )(sharedMemory + sizeOffset);
   
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;
    int cliquePartition  = warpId*pSize;
    int offsetPartition = warpId*(pSize/(k-1)+1);
    int candidatePartition = warpId*cpSize;
    int maskPartition = WarpId*cpSize*maxBitMask; 

    
    for(int i =warpId; i < totalTasks ; i+= totalWarps ){

        int start = levelData.offset[i];
        int totalCandidates = levelData.offset[i+1]- start;

        for(int iter = 0; iter <totalCandidates; iter ++){
            int candidate = levelData.candidate[start + iter];
            if(laneId==0){
                counter[threadId.x/warpSize] = 0;
            }
    
            __syncwarp();

            int candidateOffset =  candidatePartition + levelData.offsetPartition[offsetPartition + levelData.count[warpId+1]];

            
            int degree = D.degree[candidate];
            for(int j = laneId; j< degree; j+= warpSize ){
                if(label[neigh] == iterK){
                    label[neigh] = iterK-1;
                    ui loc = atomicAdd( &counter[threadId.x/warpSize], 1);
                    levelData.candidatesPartition[writeOffset + loc] = neigh;
    
                }
    
            }
            __syncwarp();
            if(laneId == 0 && counter[threadId.x/warpSize] > 0){
                levelData.partialCliquesPartition[cliquePartition + levelData.count[warpId+1] * (k-1) + level ] = vertex;
                levelData.count[warpId+1] +=1;
                levelData.offsetPartition[offsetPartition + levelData.count[warpId+1]] =
                    levelData.offsetPartition[offsetPartition + levelData.count[warpId+1] - 1] +counter[threadId.x/warpSize];
            }

        }        
    }

    __syncwarp();

    int totalTasks = levelData.count[warpId+1];

    for(int iter = 0; iter < totalTasks; iter++){
        int start = candidatePartition + levelData.offsetPartition[offsetPartition + iter];
        int end =  candidatePartition + levelData.offsetPartition[offsetPartition + iter + 1];
        int total = end-start;
        for(int i = laneId; i < total; i+=warpSize){
            int candidate = levelData.candidatesPartition[start+i]; 
            int neighOffset = D.offset[candidate];
            int degree = D.degree[candidate];

            int numBitmasks = (degree + 31) / 32;

            for (int bitmaskIndex = 0; bitmaskIndex < numBitmasks; bitmaskIndex++) {
                unsigned int bitmask = 0; // Initialize bitmask to 0

                // Iterate over the current chunk of 32 neighbors
                int startNeighbor = bitmaskIndex * 32;
                int endNeighbor = min(startNeighbor + 32, degree);
                for (int j = startNeighbor; j < endNeighbor; j++) {
                    if (label[D.neighbors[neighOffset + j]] == k - 1) {
                        bitmask |= (1 << (j - startNeighbor)); // Set the bit for valid neighbors
                    }
                }

                D.validNeighMaskPartition[ maskPartition + (levelData.offsetPartition[offsetPartition + iter] + i)*maxBitMask + bitmaskIndex] = bitmask;
            
            }


        }

    }

}

__global__ void writeFinalCliques(deviceDAGpointer D, cliqueLevelDataPointer levelData, deviceCliquesPointer trie, ui &label,ui k,ui iterK, ui n, ui m,ui psize, ui cpSize, ui maxBitMask,ui totalTasks, ui level, ui totalWarps){

    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui *counter = (ui * )(sharedMemory + sizeOffset);
   
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;
    int cliquePartition  = warpId*pSize;
    int offsetPartition = warpId*(pSize/(k-1)+1);
    int candidatePartition = warpId*cpSize;
    int maskPartition = WarpId*cpSize*maxBitMask; 

    
    for(int i =warpId; i < totalTasks ; i+= totalWarps ){

        int start = levelData.offset[i];
        int totalCandidates = levelData.offset[i+1]- start;

        for(int iter = 0; iter <totalCandidates; iter ++){
            int candidate = levelData.candidate[start + iter];
            if(laneId==0){
                counter[threadId.x/warpSize] = 0;
            }
    
            __syncwarp();

            int candidateOffset =  candidatePartition + levelData.offsetPartition[offsetPartition + levelData.count[warpId+1]];

            
            int degree = D.degree[candidate];
            for(int j = laneId; j< degree; j+= warpSize ){
                if(label[neigh] == iterK){
                    label[neigh] = iterK-1;
                    ui loc = atomicAdd( &counter[threadId.x/warpSize], 1);
                    levelData.candidatesPartition[writeOffset + loc] = neigh;
    
                }
    
            }
            __syncwarp();
            if(laneId == 0 && counter[threadId.x/warpSize] > 0){
                levelData.partialCliquesPartition[cliquePartition + levelData.count[warpId+1] * (k-1) + level ] = vertex;
                levelData.count[warpId+1] +=1;
                levelData.offsetPartition[offsetPartition + levelData.count[warpId+1]] =
                    levelData.offsetPartition[offsetPartition + levelData.count[warpId+1] - 1] +counter[threadId.x/warpSize];
            }

        }        
    }

    __syncwarp();

    int totalTasks = levelData.count[warpId+1];

    for(int iter = 0; iter < totalTasks; iter++){
        int start = candidatePartition + levelData.offsetPartition[offsetPartition + iter];
        int end =  candidatePartition + levelData.offsetPartition[offsetPartition + iter + 1];
        int total = end-start;
        for(int i = laneId; i < total; i+=warpSize){
            int candidate = levelData.candidatesPartition[start+i]; 
            int neighOffset = D.offset[candidate];
            int degree = D.degree[candidate];

            int numBitmasks = (degree + 31) / 32;

            for (int bitmaskIndex = 0; bitmaskIndex < numBitmasks; bitmaskIndex++) {
                unsigned int bitmask = 0; // Initialize bitmask to 0

                // Iterate over the current chunk of 32 neighbors
                int startNeighbor = bitmaskIndex * 32;
                int endNeighbor = min(startNeighbor + 32, degree);
                for (int j = startNeighbor; j < endNeighbor; j++) {
                    if (label[D.neighbors[neighOffset + j]] == k - 1) {
                        bitmask |= (1 << (j - startNeighbor)); // Set the bit for valid neighbors
                    }
                }

                D.validNeighMaskPartition[ maskPartition + (levelData.offsetPartition[offsetPartition + iter] + i)*maxBitMask + bitmaskIndex] = bitmask;
            
            }


        }

    }

}





