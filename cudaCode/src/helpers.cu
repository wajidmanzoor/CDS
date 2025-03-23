#include "../utils/cuda_utils.cuh"
#include "../inc/helpers.cuh"

//TODO: remove unused variables 
//TODO: Make code a bit more structured and clean
//TODO: Try to find more ways to optimize


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
        
        for (int offset = warpSize / 2; offset > 0; offset /= 2) {
            count += __shfl_down_sync(0xFFFFFFFF, count, offset);
        }

        if(laneId == 0) {
            D.degree[i] = count;
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
        }
        __syncwarp();
        ui start = G.offset[i];
        ui end = G.offset[i+1];
        ui total = end - start;
        ui neigh;
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

        int start = candidateOffset;

        for(int j = laneId; j < counter[threadIdx.x / warpSize]; j += warpSize) {
            int candidate = levelData.candidatesPartition[start + j];
            int neighOffset = D.offset[candidate];
            int degree = D.degree[candidate];

            int numBitmasks = (degree + 31) / 32;

            for (int bitmaskIndex = 0; bitmaskIndex < numBitmasks; bitmaskIndex++) {
                ui bitmask = 0; // Initialize bitmask to 0

                // Iterate over the current chunk of 32 neighbors
                int startNeighbor = bitmaskIndex * 32;
                int endNeighbor = min(startNeighbor + 32, degree);
                for (int x = startNeighbor; x < endNeighbor; x++) {


                    if (label[warpId*n + D.neighbors[neighOffset + x]] == k - 1) {
                        bitmask |= (1 << (x - startNeighbor)); // Set the bit for valid neighbors


                    }
                }



                levelData.validNeighMaskPartition[maskPartition + (levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]-1]+j) * maxBitMask + bitmaskIndex] = bitmask;
            }
        }

        __syncwarp();

        for(int i = laneId; i<n;i+=warpSize){
          label[warpId*n + i] = k;
        }

       __syncwarp();
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

    for(int iter = 0; iter < totalTasks; iter++){
        int start = candidatePartition + levelData.offsetPartition[offsetPartition + iter];
        int end = candidatePartition + levelData.offsetPartition[offsetPartition + iter+ 1];
        int total = end-start;

        int writeOffset = levelData.temp[warpId] + levelData.offsetPartition[offsetPartition + iter];
        for(int i = laneId; i < total; i+=warpSize){
           ui candidate = levelData.candidatesPartition[start + i];
            levelData.candidates[writeOffset+ i] = levelData.candidatesPartition[start + i];
            
            int totalMasks = (D.degree[candidate]+31)/32;
            for(int j =0; j < totalMasks; j++){
                levelData.validNeighMask[(writeOffset+i)*maxBitMask + j ] =
                levelData.validNeighMaskPartition[maskPartition + (levelData.offsetPartition[offsetPartition + iter] + i)*maxBitMask + j];
            }

        }

        if(laneId< level+1 ){

                levelData.partialCliques[levelData.count[warpId]*(k-1)+ iter*(k-1) + laneId] = levelData.partialCliquesPartition[cliquePartition + iter * (k-1) + laneId];
          }

        __syncwarp();

        if(laneId==0){

            levelData.offset[levelData.count[warpId] + iter + 1] = levelData.temp[warpId]+levelData.offsetPartition[offsetPartition + iter+ 1];

        }
        __syncwarp();


    }

}


__global__ void listMidCliques(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui *label, ui k,ui iterK, ui n, ui m,ui pSize, ui cpSize, ui maxBitMask,ui totalTasks, ui level, ui totalWarps){

    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui *counter = (ui * )(sharedMemory + sizeOffset);

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;
    int cliquePartition  = warpId*pSize;
    int offsetPartition = warpId*(pSize/(k-1)+1);
    int candidatePartition = warpId*cpSize;
    int maskPartition = warpId*cpSize*maxBitMask;


    for(int i =warpId; i < totalTasks ; i+= totalWarps ){

        int start = levelData.offset[i];
        int totalCandidates = levelData.offset[i+1]- start;

        for(int iter = 0; iter <totalCandidates; iter ++){
            int candidate = levelData.candidates[start + iter];
            if(laneId==0){
                counter[threadIdx.x/warpSize] = 0;
            }
            __syncwarp();

            int degree = D.degree[candidate];
            int neighOffset = D.offset[candidate];

            int writeOffset = candidatePartition + levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]];
            for(int j = laneId; j< degree; j+= warpSize ){
                int iterBitMask = j/warpSize;
                int bitPos = j%32;
                int neighBitMask = levelData.validNeighMask[start*maxBitMask + iter + iterBitMask];
                ui neigh = D.neighbors[neighOffset + j];
                if( (label[warpId*n + neigh] == iterK) && (neighBitMask & (1 << bitPos )) ){
                    label[warpId*n + neigh] = iterK-1;
                    ui loc = atomicAdd( &counter[threadIdx.x/warpSize], 1);

                    levelData.candidatesPartition[writeOffset + loc] = neigh;

                }

            }
            __syncwarp();
            if(laneId == 0 && counter[threadIdx.x/warpSize] > 0){
                levelData.partialCliquesPartition[cliquePartition + levelData.count[warpId+1] * (k-1) + level ] = candidate;
                for(int l =0; l<level; l++){
                  levelData.partialCliquesPartition[cliquePartition + levelData.count[warpId+1] * (k-1) + l ] = levelData.partialCliques[i*(k-1)+l];
                }
                levelData.count[warpId+1] +=1;
                levelData.offsetPartition[offsetPartition + levelData.count[warpId+1]] =
                    levelData.offsetPartition[offsetPartition + levelData.count[warpId+1] - 1] +counter[threadIdx.x/warpSize];
            }

            __syncwarp();
            int start = writeOffset;

            for(int j = laneId; j < counter[threadIdx.x / warpSize]; j += warpSize) {
                int cand = levelData.candidatesPartition[start + j];
                int neighOffset = D.offset[cand];
                int degree = D.degree[cand];
    
                int numBitmasks = (degree + 31) / 32;
    
                for (int bitmaskIndex = 0; bitmaskIndex < numBitmasks; bitmaskIndex++) {
                    ui bitmask = 0; // Initialize bitmask to 0
    
                    // Iterate over the current chunk of 32 neighbors
                    int startNeighbor = bitmaskIndex * 32;
                    int endNeighbor = min(startNeighbor + 32, degree);
                    for (int x = startNeighbor; x < endNeighbor; x++) {
    
    
                        if (label[warpId*n + D.neighbors[neighOffset + x]] == iterK - 1) {
                            bitmask |= (1 << (x - startNeighbor)); // Set the bit for valid neighbors
    
    
                        }
                    }
    
    
    
                    levelData.validNeighMaskPartition[maskPartition + (levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]-1]+j) * maxBitMask + bitmaskIndex] = bitmask;
                }
            }
    
            __syncwarp();
    
            for(int i = laneId; i<n;i+=warpSize){
              label[warpId*n + i] = iterK;
            }
    
           __syncwarp();
        }

    }
    
}


__global__ void writeFinalCliques(deviceGraphPointers G, deviceDAGpointer D, cliqueLevelDataPointer levelData, deviceCliquesPointer cliqueData, ui *globalCounter,ui k,ui iterK, ui n, ui m,ui pSize, ui cpSize, ui maxBitMask,ui trieSize,ui totalTasks, ui level, ui totalWarps){
    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui *counter = (ui * )(sharedMemory + sizeOffset);
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    
    for(int i =warpId; i < totalTasks ; i+= totalWarps ){

        int start = levelData.offset[i];
        int totalCandidates = levelData.offset[i+1]- start;

        for(int iter = 0; iter <totalCandidates; iter ++){
            int candidate = levelData.candidates[start + iter];
            if(laneId==0){
                counter[warpId]=0;
            }
            __syncwarp();
            int degree = D.degree[candidate];
            int neighOffset = D.offset[candidate];

            
            for(int j = laneId; j< degree; j+= warpSize ){
                int iterBitMask = j/warpSize;
                int bitPos = j%32;
                int neighBitMask = levelData.validNeighMask[start*maxBitMask + iter + iterBitMask];
                if(neighBitMask & (1 << bitPos )){

                    ui neigh = D.neighbors[neighOffset + j];
                  
                    ui loc = atomicAdd(globalCounter,1);
                    for(int ind =0; ind < k-2; ind++){
                        cliqueData.trie[trieSize * ind + loc] = levelData.partialCliques[(i)*(k-1) + ind];
                        
                    }
                    atomicAdd(&counter[warpId],1);
                    cliqueData.trie[trieSize * (k-2) + loc]  = candidate;
                    cliqueData.trie[trieSize * (k-1) + loc] = neigh;
                    cliqueData.status[loc]=1;
                    atomicAdd(&G.cliqueDegree[neigh],1);
                    atomicAdd(&G.cliqueDegree[candidate],1);

                }
                
    
            }
            __syncwarp();

            for(int j = laneId; j< k-2 ; j+= warpSize ){
                int pClique = levelData.partialCliques[i*(k-1) + j];
                atomicAdd(&G.cliqueDegree[pClique],counter[warpId]);
            }

        }        
    }

}


__global__ void sortTrieData(deviceGraphPointers G, deviceCliquesPointer cliqueData, ui totalCliques, ui t, ui k, ui totalThreads){
    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui *elements = (ui * )(sharedMemory + sizeOffset);
    sizeOffset = k*WARPS_EACH_BLK*sizeof(ui);
    ui *degree = (ui * )(sharedMemory + sizeOffset);


    int idx = blockIdx.x * blockDim.x + threadIdx.x;


    for(int i = idx; i <totalCliques; i+=totalThreads ){

        for(int j=0;j<k;j++){
            
            elements[j] = cliqueData.trie[j*t+i];
            degree[j] = G.cliqueDegree[elements[j]];
        }

        __syncwarp();

        // Use insertion sort, as it is best for small arrays 

        for(int j=1;j<k;j++){
            ui current_element = elements[j];
            ui current_degree = degree[j];
            int ind = j-1;

            while(ind >= 0 && degree[ind] > current_degree){
                elements[ind + 1] = elements[ind];
                degree[ind + 1] = degree[ind];
                ind--;
            }
            elements[ind + 1] = current_element;
            degree[ind + 1] = current_degree;

        }

        for(int j=0;j<k;j++){
            cliqueData.trie[j*t+i] = elements[j];
    
        }

    }

}


__global__ void selectNodes(deviceGraphPointers G, ui *bufTails,ui *glBuffers, ui glBufferSize, ui n, ui level){
    __shared__ ui *glBuffer;
    __shared__ ui bufTail;

    if(threadIdx.x == 0){
        bufTail = 0;
        glBuffer = glBuffers + blockIdx.x*glBufferSize;
    }
    __syncthreads();

    ui idx = blockIdx.x * blockDim.x + threadIdx.x;
    for(ui i = idx ;i<n; i+=BLK_DIM){
      ui v = i;

      if(G.cliqueCore[v] == level){
        ui loc = atomicAdd(&bufTail, 1);
        glBuffer[loc] = v;
        }
    }
    __syncthreads();

    if(threadIdx.x == 0)
    {
        bufTails [blockIdx.x] = bufTail;
    }

}

__global__ void processNodesByWarp(deviceGraphPointers G,deviceCliquesPointer cliqueData, ui *bufTails,ui *glBuffers, ui *globalCount, ui glBufferSize, ui n, ui level, ui k, ui t){
    __shared__ ui bufTail;
    __shared__ ui *glBuffer;
    __shared__ ui base;
    ui warpId = threadIdx.x / 32;
    ui laneId = threadIdx.x % 32;
    ui regTail;
    ui i;
    if(threadIdx.x==0){
    bufTail = bufTails[blockIdx.x];
    base = 0;
    glBuffer = glBuffers + blockIdx.x*glBufferSize;
    assert(glBuffer!=NULL);
    }

    while(true){
    __syncthreads();
    if(base == bufTail) break; // all the threads will evaluate to true at same iteration
    i = base + warpId;
    regTail = bufTail;
    __syncthreads();

    if(i >= regTail) continue; // this warp won't have to do anything

    if(threadIdx.x == 0){
    base += WARPS_EACH_BLK;
    if(regTail < base )
    base = regTail;
    }
    //bufTail is incremented in the code below:
    ui v = glBuffer[i];


   __syncthreads();
    for(ui j =laneId; j<t; j+=warpSize){

        if( (v == cliqueData.trie[j]) && (cliqueData.status[j] ==1)){

            for(ui x =1;x<k;x++){
                ui u = cliqueData.trie[x*t+i];
                int a = atomicSub(&G.cliqueCore[u],1);
                if(a == level+1){
                    ui loc = atomicAdd(&bufTail, 1);
                    glBuffer[loc] = u;

                }
                if(a <= level){
                    atomicAdd(&G.cliqueCore[u], 1);
                }
            }
            cliqueData.status[i] = 0;


        }
    }

__syncthreads();
    if(threadIdx.x == 0 && bufTail>0){
    atomicAdd(globalCount, bufTail); // atomic since contention among blocks
    }
}
}

__global__ void processNodesByBlock(deviceGraphPointers G,deviceCliquesPointer cliqueData, ui *bufTails,ui *glBuffers, ui *globalCount, ui glBufferSize, ui n, ui level, ui k, ui t){
    __shared__ ui bufTail;
    __shared__ ui *glBuffer;
    __shared__ ui base;

    ui regTail;
    ui i;
    if(threadIdx.x==0){
    bufTail = bufTails[blockIdx.x];
    base = 0;
    glBuffer = glBuffers + blockIdx.x*glBufferSize;
    assert(glBuffer!=NULL);
    }

    while(true){
    __syncthreads();
    if(base == bufTail) break; // all the threads will evaluate to true at same iteration
    i = base + blockIdx.x;
    regTail = bufTail;
    __syncthreads();

    if(i >= regTail) continue; // this warp won't have to do anything

    if(threadIdx.x == 0){
    base += 1;
    if(regTail < base )
    base = regTail;
    }
    //bufTail is incremented in the code below:
    ui v = glBuffer[i];


   __syncthreads();
    int idx = threadIdx.x ;
    for(ui j = idx; j<t; j+= BLK_DIM){
        if( (v = cliqueData.trie[j]) && (cliqueData.status[j]==1)){
            for(ui x =1;x<k;x++){
                ui u = cliqueData.trie[x*t+i];
                int a = atomicSub(&G.cliqueCore[u], 1);
                if(a == level+1){
                    ui loc = atomicAdd(&bufTail, 1);
                    glBuffer[loc] = u;

                }
                if(a <= level){
                    atomicAdd(&G.cliqueCore[u], 1);
                }
            }
            cliqueData.status[i] = 0;


        }

    }


    __syncthreads();

    if(threadIdx.x == 0 && bufTail>0){
    atomicAdd(globalCount, bufTail); // atomic since contention among blocks
    }
}
}


__global__ void locateDensetCore(){
    
}