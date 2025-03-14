#include "../utils/cuda_utils.cuh"
#include "../inc/helpers.cuh"



__global__ void generateDegreeDAG(deviceGraphPointers G, deviceDAGpointer D, ui *listingOrder, ui n, ui m, ui totalWarps) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(ui i = warpId; i < n; i += totalWarps) {
        ui start = G.offset[i];
        ui end = G.offset[i+1];
        ui total = end - start;
        ui neigh;
        int count = 0;
        for(int j = laneId; j < total; j += warpSize) {
            neigh = G.neighbors[start + j];
            if(listingOrder[i] < listingOrder[neigh]) {
                count++;
            }
        }
        unsigned mask = __ballot_sync(0xFFFFFFFF, count > 0);
        int total_count = __popc(mask);

        if(laneId == 0) {
            D.degree[i] = total_count;
        }
    }
}

__global__ void generateNeighborDAG(deviceGraphPointers G, deviceDAGpointer D, ui *listingOrder, ui n, ui m, ui totalWarps) {

     extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui *counter = (ui *)(sharedMemory + sizeOffset);

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;
   
    for(ui i = warpId; i < n; i += totalWarps) {
        if(laneId==0){
          counter[threadIdx.x / warpSize] = D.offset[i];
          //printf("warp %d  counter[%d] = %d\n", i, threadIdx.x / warpSize, counter[threadIdx.x / warpSize]);
        }
        __syncwarp();
        ui start = G.offset[i];
        ui end = G.offset[i+1];
        ui total = end - start;
        ui neigh;
        //printf("warp %d  start %d end %d total %d\n", i, start, end, total);
        for(int j = laneId; j < total; j += warpSize) {
            neigh = G.neighbors[start + j];

            if(listingOrder[i] < listingOrder[neigh]) {
                int loc = atomicAdd(&counter[threadIdx.x / warpSize], 1);
                D.neighbors[loc] = neigh;

            }
        }
      __syncwarp();
    }
}

__global__ void listIntialCliques(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui *label, ui k, ui n, ui m, ui psize, ui cpSize, ui maxBitMask, ui level, ui totalWarps) {
    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui *counter = (ui *)(sharedMemory + sizeOffset);

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;
    int cliquePartition = warpId * psize;
    int offsetPartition = warpId * (psize / (k-1) + 1);
    int candidatePartition = warpId * cpSize;
    int maskPartition = warpId * cpSize * maxBitMask;

    for(int i = warpId; i < n; i += totalWarps) {

        ui vertex = i;
        ui neighOffset = D.offset[vertex];
        if(laneId == 0) {
            counter[threadIdx.x / warpSize] = 0;
        }

        __syncwarp();
        int candidateOffset = candidatePartition + levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]];
        

        for(int j = laneId; j < D.degree[vertex]; j += warpSize) {
            ui neigh = D.neighbors[neighOffset + j];
            if(label[warpId*n + neigh] == k) {
                label[warpId*n + neigh] = k - 1;
                ui loc = atomicAdd(&counter[threadIdx.x / warpSize], 1);
                levelData.candidatesPartition[candidateOffset + loc] = neigh;
            }
        }
        __syncwarp();
        if(laneId == 0 && counter[threadIdx.x / warpSize] > 0) {
            levelData.partialCliquesPartition[cliquePartition + levelData.count[warpId + 1] * (k-1) + level] = vertex;
            levelData.count[warpId + 1] += 1;
            levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]] =
                levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1] - 1] + counter[threadIdx.x / warpSize];
        }
         __syncwarp();
        for(int i = laneId; i<n;i+=32){
          label[warpId*n + i] = k;
        }
    }

    __syncwarp();

    int totalTasks = levelData.count[warpId + 1];

    for(int iter = 0; iter < totalTasks; iter++) {
        int start = candidatePartition + levelData.offsetPartition[offsetPartition + iter];
        int end = candidatePartition + levelData.offsetPartition[offsetPartition + iter + 1];
        int total = end - start;
        for(int i = laneId; i < total; i += warpSize) {
            int candidate = levelData.candidatesPartition[start + i];
            int neighOffset = D.offset[candidate];
            int degree = D.degree[candidate];

            int numBitmasks = (degree + 31) / 32;

            for (int bitmaskIndex = 0; bitmaskIndex < numBitmasks; bitmaskIndex++) {
                unsigned int bitmask = 0; // Initialize bitmask to 0

                // Iterate over the current chunk of 32 neighbors
                int startNeighbor = bitmaskIndex * 32;
                int endNeighbor = min(startNeighbor + 32, degree);
                for (int j = startNeighbor; j < endNeighbor; j++) {
                    if (label[warpId*n + D.neighbors[neighOffset + j]] == k - 1) {
                        bitmask |= (1 << (j - startNeighbor)); // Set the bit for valid neighbors
                    }
                }

                levelData.validNeighMaskPartition[maskPartition + (levelData.offsetPartition[offsetPartition + iter] + i) * maxBitMask + bitmaskIndex] = bitmask;
            }
        }
    }
   
}


__global__ void flushParitions(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui pSize, ui cpSize, ui k, ui maxBitMask, ui level, ui totalWarps){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;
    int cliquePartition = warpId * pSize;
    int offsetPartition = warpId * (pSize / (k-1) + 1);
    int candidatePartition = warpId * cpSize;
    int maskPartition = warpId * cpSize * maxBitMask;

    int totalTasks = levelData.count[warpId+1] - levelData.count[warpId];
    //printf("warp %d totalTasks %d, offsetPartition %d offset part %d \n", warpId, totalTasks,candidatePartition,offsetPartition);

    for(int iter = 0; iter < totalTasks; iter++){
        int start = candidatePartition + levelData.offsetPartition[offsetPartition + iter];
        int end = candidatePartition + levelData.offsetPartition[offsetPartition + iter+ 1];
        int total = end-start;
        //printf("warp %d start %d end %d total %d\n", warpId, start, end, total);
        int writeOffset = levelData.temp[warpId] + levelData.offsetPartition[offsetPartition + iter];
        for(int i = laneId; i < total; i+=warpSize){
           ui candidate = levelData.candidatesPartition[start + i];
            levelData.candidates[writeOffset+ i] = levelData.candidatesPartition[start + i];
            if(i< level){
          
                levelData.partialCliques[levelData.count[warpId]*(k-1)+ iter*(k-1) + i] = levelData.partialCliquesPartition[cliquePartition + iter * (k-1) + i];
                printf("warp Id %d lane id %d count %d iter % d pclique %d level %d write loc %d \n",warpId,i,levelData.count[warpId],iter,levelData.partialCliquesPartition[cliquePartition + iter * (k-1) + i],i,levelData.count[warpId]+ iter*(k-1) + i);
            }
            int totalMasks = (D.degree[candidate]+31)/32;
            for(int j =0; j < totalMasks; j++){
                levelData.validNeighMask[(writeOffset+i)*maxBitMask + j ] = 
                levelData.validNeighMaskPartition[maskPartition + (levelData.offsetPartition[offsetPartition + iter] + i)*maxBitMask + j];
            }

        }

        __syncwarp();

        if(laneId==0){
            levelData.offset[levelData.count[warpId] + iter] = levelData.temp[warpId]+levelData.offsetPartition[offsetPartition + iter+ 1];

        }


    }

}



__global__ void listMidCliques(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui *label, ui k,ui iterK, ui n, ui m,ui psize, ui cpSize, ui maxBitMask,ui totalTasks, ui level, ui totalWarps){

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
                ui neigh = D.neighbors[neighOffset + j];

                if(label[warpId*n + neigh] == iterK){
                    label[warpId*n + neigh] = iterK-1;
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
                    if (label[warpId*n + D.neighbors[neighOffset + j]] == k - 1) {
                        bitmask |= (1 << (j - startNeighbor)); // Set the bit for valid neighbors
                    }
                }

                D.validNeighMaskPartition[ maskPartition + (levelData.offsetPartition[offsetPartition + iter] + i)*maxBitMask + bitmaskIndex] = bitmask;
            
            }


        }

    }

}

__global__ void writeFinalCliques(deviceGraphPointers G, deviceDAGpointer D, cliqueLevelDataPointer levelData, deviceCliquesPointer cliqueData, ui *globalCounter,ui k,ui iterK, ui n, ui m,ui psize, ui cpSize, ui maxBitMask,ui totalTasks, ui level, ui totalWarps){
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
            int candidateOffset =  candidatePartition + levelData.offsetPartition[offsetPartition + levelData.count[warpId+1]];

            
            int degree = D.degree[candidate];
            int iterMask = 0;
            if(laneId==0){
                counter[warpId]=0;
            }
            for(int j = laneId; j< degree; j+= warpSize ){

                int chunkEnd = min( (iterMaks + 1)*32, degree);
                ui chunkMask = levelData.validNeighMask[(candidateOffset + i)*maxBitMask + iterMask];
                ui lshift = j%32;
                if(chunkMask & (1 << lshift)){
                    ui neigh = D.neighbors[neighOffset + j];
                    
                    ui loc = atomicAdd(globalCounter,1);
                    atomicAdd(&counter[warpId],1);
                    for(int ind =0; ind < k-2; ind++){
                        cliqueData.trie[trieSize * ind + loc] = levelData.partialCliques[levelData.count[warpId]+ iter*(k-1) + ind];
                    }
                    cliqueData.trie[trieSize * (k-2) + loc]  = candidate;
                    cliqueData.trie[trieSize * (k-1) + loc] = neigh;
                    cliqueData.status[loc]=1;
                    atomicAdd(&G.cliqueDegree[neigh],1);
                    atomicAdd(&G.cliqueDegree[candidate],1);

                }
                
    
            }
            __syncwarp();

            for(int j = laneId; j< k-2 ; j+= warpSize ){
                int pClique = levelData.partialCliques[levelData.count[warpId]+ iter*(k-1) + j];
                atomicAdd(&G.cliqueDegree[pClique],counter[warpId]);
            }

        }        
    }

}




__global__ void sortTrieData(deviceGraphPointers G, deviceCliquesPointer cliqueData, ui t, ui k, ui totalThreads){

    int idx = blockIdx.x * blockDim.x + threadIdx.x;


    for(int i = idx; i <t; i+=totalThreads ){
        ui elements[k];
        ui degree[k];
        for(int j=0;j<k;j++){
            elements[j] = cliqueData.trie[j*t+i];
            degree[j] = G.cliqueDegree[element];
        }

        // Use insertion sort, as it is best for small arrays 

        for(j=1;j<k;j++){
            ui insert_index = j;
            int current_element = elements[j];
            int current_degree = degree[current_element];
            ui current = degree[element[j]];
            for(int ind = j-1; ind>=0; ind-- ){
                if(degree[elements[ind]] > ){current_degree
                    elements[ind+1] = elements[ind];
                    insert_index = ind;

                }else{
                    break;
                }

            }
            elements[insert_index] = current_element;

        }

        for(int j=0;j<k;j++){
            cliqueData.trie[j*t+i] = elements[j];
    
        }

    }

}
