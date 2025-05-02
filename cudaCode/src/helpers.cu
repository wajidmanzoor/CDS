#include "../utils/cuda_utils.cuh"
#include "../inc/helpers.cuh"

//TODO: remove unused variables 
//TODO: Make code a bit more structured and clean
//TODO: Try to find more ways to optimize

__device__ double fact(ui k){
    double res = 1;
    int i = k;
    while(i>1){
        res= res*i;
        i--;

    }
    return res;
}

__device__ double power(ui totalCliques, double p) {
    return powf((double)totalCliques, p );
}

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
                    cliqueData.status[loc]=-1;
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

__global__ void processNodesByWarp(deviceGraphPointers G,deviceCliquesPointer cliqueData, ui *bufTails,ui *glBuffers, ui *globalCount, ui glBufferSize, ui n, ui level, ui k, ui t, ui tt){
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


   __syncwarp();
    for(ui j =laneId; j<tt; j+=warpSize){
    //printf("warpId %d laneId %u vertex %u check %d t %d \n",warpId,i,v,cliqueData.trie[j],t);
        if(cliqueData.status[j] == -1){
            bool found = false;
            ui w =0;
            while(w<k){
                if(cliqueData.trie[w*t+j] == v){
                found = true;
                break;
                }
                w++;
            }
            if(found){

                for(ui x =0;x<k;x++){

                    if(x==w) continue;
                    ui u = cliqueData.trie[x*t+j];
                    int a = atomicSub(&G.cliqueCore[u],1);
                    if(a == level+1){
                        ui loc = atomicAdd(&bufTail, 1);
                        glBuffer[loc] = u;

                    }
                    if(a <= level){
                        atomicAdd(&G.cliqueCore[u], 1);
                    }

                }
                cliqueData.status[j] = level;


            }
        }
    }

    __syncwarp();
    if(laneId == 0 && bufTail>0){
      atomicAdd(globalCount, 1); // atomic since contention among blocks
    }
}
}

__global__ void processNodesByBlock(deviceGraphPointers G,deviceCliquesPointer cliqueData, ui *bufTails,ui *glBuffers, ui *globalCount, ui glBufferSize, ui n, ui level, ui k, ui t, ui tt){
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
        ui idx = threadIdx.x;

        for(ui j = idx; j<tt; j+= BLK_DIM){
            
            if(cliqueData.status[j]==-1){

                bool found = false;
                ui w =0;
                while(w<k){
                    if(cliqueData.trie[w*t+j] == v){
                    found = true;
                    break;
                    }
                    w++;
                }

                if(found){
                    for(ui x =0;x<k;x++){
                        if(x==w) continue;

                        ui u = cliqueData.trie[x*t+j];
                        int a = atomicSub(&G.cliqueCore[u], 1);
                        if(a == level+1){
                            ui loc = atomicAdd(&bufTail, 1);
                            glBuffer[loc] = u;

                        }
                        if(a <= level){
                            atomicAdd(&G.cliqueCore[u], 1);
                        }
                    }
                    cliqueData.status[j] = level;


                }
            }
        }


        __syncthreads();

        if(threadIdx.x == 0 && bufTail>0){
            atomicAdd(globalCount, 1); // atomic since contention among blocks
        }
    }
}

__global__ void generateDensestCore(deviceGraphPointers G, densestCorePointer densestCore,ui *globalCount, ui n, ui maxDensityCore, ui totalWarps){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(ui i = warpId; i < n; i += totalWarps){
        if(G.cliqueCore[i]>= maxDensityCore){
            ui loc;
            if(laneId==0){
                loc = atomicAdd(globalCount,1);
                densestCore.mapping[loc] = i;
            }
            loc = __shfl_sync(0xFFFFFFFF, loc, 0, 32);
            ui start = G.offset[i];
            ui end = G.offset[i+1];
            ui total = end - start;
            ui neigh;
            int count = 0;
            for(int j = laneId; j < total; j += warpSize) {
                neigh = G.neighbors[start + j];
                if(G.cliqueCore[neigh] >= maxDensityCore) {
                    count++;
                }
            }
            for (int offset = warpSize / 2; offset > 0; offset /= 2) {
                count += __shfl_down_sync(0xFFFFFFFF, count, offset);
            }
            if(laneId == 0) {
                densestCore.offset[loc+1] = count;
            }

        }
    }
}

__global__ void generateNeighborDensestCore(deviceGraphPointers G, densestCorePointer densestCore, ui maxDensityCore, ui totalWarps) {

    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui *counter = (ui *)(sharedMemory + sizeOffset);

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(ui i = warpId; i < (*densestCore.n); i += totalWarps) {
        if(laneId==0){
          counter[threadIdx.x / warpSize] = densestCore.offset[i];
        }
        __syncwarp();
        ui vertex = densestCore.mapping[i];
        ui start = G.offset[vertex];
        ui end = G.offset[vertex+1];
        ui total = end - start;
        ui neigh;
        for(int j = laneId; j < total; j += warpSize) {
            neigh = G.neighbors[start + j];

            if(G.cliqueCore[neigh] >= maxDensityCore) {
                int loc = atomicAdd(&counter[threadIdx.x / warpSize], 1);

                densestCore.neighbors[loc] = densestCore.reverseMap[neigh];

            }
        }
      __syncwarp();
    }
}

__global__ void pruneEdges(densestCorePointer densestCore, deviceCliquesPointer cliqueData, ui *pruneStatus,ui t, ui tt, ui k, ui level ){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(ui i = warpId; i<tt; i+=TOTAL_WARPS){

        if(cliqueData.status[i] >= level){

            for(ui iter =0; iter< k ; iter ++){
                // v should be mapped
                ui u_ = ((iter)%k)*t+i;
                ui u  =  densestCore.reverseMap[cliqueData.trie[u_]];
                for(ui j = 0; j < k; j++){
                    ui v_ = ((j)%k)*t+i;


                    if(v_!=u_){
                        int v = densestCore.reverseMap[cliqueData.trie[v_]];


                        // Update u-v edge status
                        ui start = densestCore.offset[u];
                        ui end = densestCore.offset[u+1];
                        ui total = end-start;

                        for(ui ind = laneId; ind < total; ind +=WARPSIZE){
                            int neigh = densestCore.neighbors[start+ind];
                            if(neigh == v){
                              atomicCAS(&pruneStatus[start+ind], 1, 0);
                            }

                        }

                    }
                }




        }

        }
    }
}

__global__ void generateDegreeAfterPrune(densestCorePointer densestCore ,ui *pruneStatus, ui *newOffset, ui n, ui m, ui totalWarps) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(ui i = warpId; i < n; i += totalWarps) {
        ui start = densestCore.offset[i];
        ui end = densestCore.offset[i+1];
        ui total = end - start;
        int count = 0;
        for(int j = laneId; j < total; j += warpSize) {
            if(!pruneStatus[start + j]) {
                count++;
            }
        }

        for (int offset = warpSize / 2; offset > 0; offset /= 2) {
            count += __shfl_down_sync(0xFFFFFFFF, count, offset);
        }

        if(laneId == 0) {
            newOffset[i+1] = count;
        }
    }
}

__global__ void generateNeighborAfterPrune(densestCorePointer densestCore ,ui *pruneStatus, ui *newOffset, ui *newNeighbors,ui n, ui m, ui totalWarps) {

    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui *counter = (ui *)(sharedMemory + sizeOffset);

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(ui i = warpId; i < n; i += totalWarps) {
        if(laneId==0){
          counter[threadIdx.x / warpSize] = newOffset[i];
        }
        __syncwarp();
        ui start = densestCore.offset[i];
        ui end = densestCore.offset[i+1];
        ui total = end - start;
        ui neigh;
        for(int j = laneId; j < total; j += warpSize) {
            neigh = densestCore.neighbors[start + j];

            if(!pruneStatus[start + j]) {
                int loc = atomicAdd(&counter[threadIdx.x / warpSize], 1);
                newNeighbors[loc] = neigh;

            }
        }
      __syncwarp();
    }
}

__global__ void componentDecomposek(deviceComponentPointers conComp, devicePrunedNeighbors prunedNeighbors, ui *changed, ui n, ui m, ui totalWarps) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;
    bool threadChanged = false;

    for(ui i = warpId; i < n; i += totalWarps) {
        ui currentComp = conComp.components[i];
        ui start = prunedNeighbors.newOffset[i];
        ui end = prunedNeighbors.newOffset[i+1];
        ui total = end - start;
        //printf("warpid %d laneId %d start %d end %d total %d cc %d \n",warpId,laneId,start,end,total,currentComp);

        ui minNeighComp = currentComp;

        for (ui j = laneId; j < total; j += warpSize) {
            ui neighComp = conComp.components[prunedNeighbors.newNeighbors[start+j]];
            minNeighComp = min(minNeighComp, neighComp);
            //printf("warp Id %d laneid %d cc %d nc %d mc %d \n",warpId,laneId,currentComp,neighComp,minNeighComp);
        }

        for (int offset = warpSize/2; offset > 0; offset /= 2) {
            ui temp = __shfl_down_sync(0xFFFFFFFF, minNeighComp, offset);
            minNeighComp = min(minNeighComp, temp);
        }

        if (laneId == 0) {
            if ( minNeighComp < currentComp) {
                conComp.components[i] = minNeighComp;
                threadChanged = true;
            }
        }

        __syncwarp();
    }

    bool warpChanged = __any_sync(0xFFFFFFFF, threadChanged);
    if (warpChanged && laneId == 0) {
        atomicAdd(changed, 1);
    }
}

__global__ void getConnectedComponentStatus(deviceComponentPointers conComp,deviceCliquesPointer cliqueData, densestCorePointer densestCore, ui *compCounter, ui t, ui tt, ui k,ui maxCore, ui totalThreads){

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    for(ui i =idx; i<tt;i +=totalThreads){
        if(cliqueData.status[i]>=maxCore){
          int comp = INT_MAX;

          for(ui x=0;x<k;x++){
            ui vertex =densestCore.reverseMap[cliqueData.trie[x*t + i]];
            comp = min(comp,conComp.components[vertex]);

            cliqueData.trie[x*t + i] = vertex;


          }
          cliqueData.status[i] = comp;
          atomicAdd(&compCounter[comp+1],1);

        
      }else{
        cliqueData.status[i] = -1;
      }
    }

}

__global__ void rearrangeCliqueData(deviceComponentPointers conComp,deviceCliquesPointer cliqueData, deviceCliquesPointer finalCliqueData,densestCorePointer densestCore, ui *compCounter,ui *counter,ui t, ui tt, ui k,ui totaLCliques, ui totalThreads){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    for(ui i =idx; i<tt;i +=totalThreads){
        
        int comp =  cliqueData.status[i];

        if(comp>-1){
          ui loc; 
          for(ui j=0;j<k;j++){
              ui vertex = cliqueData.trie[j*t + i];
              ui offset = compCounter[comp];
              if(j==0){
                  loc = atomicAdd(&counter[comp],1);
                  finalCliqueData.status[offset + loc ] = comp;


              }
              finalCliqueData.trie[offset + j*totaLCliques + loc ] = vertex;
          }
        }
        
    }

}
__global__ void getLbUbandSize(deviceComponentPointers conComp, ui *compCounter, ui *lowerBound, ui *upperBound, ui *ccOffset,  ui *neighborSize, ui totalComponenets, ui k, ui maxDensity){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx==0){
        ccOffset[idx]=0;
    }

    for(ui i = idx; i<totalComponenets; i+=TOTAL_THREAD ){
        ui totalCliques = compCounter[i+1] - compCounter[i];
        ui totalSize = conComp.componentOffset[i+1] -  conComp.componentOffset[i];
        double lb = (double) (totalCliques)/totalSize;
        atomicMax(lowerBound, lb);

        double den = power(fact(k),1.0/k);
        double num = power(totalCliques, (k-1.0)/k);
        double up = max(maxDensity, num/dem);

        upperBound[i] = up;

        if(up>lb){
            ccOffset[i+1] = totalCliques + totalSize + 2;
            atomicAdd(neighborSize, 2*totalCliques*k + 4*totalSize);


        }
        else{
            ccOffset[i+1] = 0;
        }
    }


}

__global__ void createFlowNetwork(deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp, densestCorePointer densestCore, deviceCliquesPointer finalCliqueData, ui *compCounter,ui *counter,ui *upperBound, ui *offset,ui *rank, ui n, ui m, ui totalWarps, int totalComponents, ui k, ui lb) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(ui i = warpId; i < totalComponents; i += totalWarps){

        if(upperBound[i]>lb){
            ui alpha = (upperBound[i]+lb)/2;
            ui offsetLoc;
            ui start = conComp.componentOffset[i];
            ui end = conComp.componentOffset[i+1];
            ui total = end - start;
            ui cliqueStart = compCounter[i];
            ui cliqueEnd = compCounter[i+1];
            ui totalCliques = cliqueEnd-cliqueStart;

            ui offset = offset[rank[i]];

            for (ui j = laneId; j < totalCliques; j += warpSize){
                for(ui x =0; x <k; x++){
                    ui vertex = cliqueData.finalCliqueData[cliqueStart + x*t+j];
                    flowNetwork.toEdge[offset + j* k +x ] = vertex;
                    flowNetwork.capacity[offset + j* k +x ] = k-1;
                    flowNetwork.flow[offset + j* k +x ] = k-1;

                    ui loc = atomicAdd(&counter[vertex],1);
                    flowNetwork.toEdge[offset + totalCliques + vertex + loc] = j;
                    flowNetwork.capacity[offset + totalCliques + vertex + loc] = 1;
                    flowNetwork.flow[offset + totalCliques + vertex + loc ] = 1;
                }
            }

            for (ui j = laneId; j < total; j += warpSize){
                ui vertex = conComp.mapping[offset+j];
                flowNetwork.toEdge[offset + totalCliques + total + 1] = vertex;
                flowNetwork.capacity[offset + totalCliques + total + 1] = densestCore.cliqueDegree[vertex];
                flowNetwork.flow[offset + totalCliques + total + 1] = densestCore.cliqueDegree[vertex];

                ui loc = atomicAdd(&counter[vertex],1);
                flowNetwork.toEdge[offset + totalCliques + vertex + loc] = total+totalCliques;
                flowNetwork.capacity[offset + totalCliques + vertex + loc] = 0;
                flowNetwork.flow[offset + totalCliques + vertex + loc ] = 0;

                flowNetwork.toEdge[offset + totalCliques + total + 2] = vertex;
                flowNetwork.capacity[offset + totalCliques + total + 1] = 0;
                flowNetwork.flow[offset + totalCliques + total + 1] = 0;

                loc = atomicAdd(&counter[vertex],1);
                flowNetwork.toEdge[offset + totalCliques + vertex + loc] = total+totalCliques+1;
                flowNetwork.capacity[offset + totalCliques + vertex + loc] = alpha*k;
                flowNetwork.flow[offset + totalCliques + vertex + loc ] = alpha*k;

            }

            for(ui j = laneId; j < totalCliques + total + 2; j += warpSize){
                if(j ==0){
                    flowNetwork.offset[need + j ] =0;
                } else if(j<totalCliques+1){
                    flowNetwork.offset[need + j ]  = j*k; 

                } else if(j<totalCliques+total+1){
                    flowNetwork.offset[need + j ] = totalCliques*k+

                }

                //Need to implement a way to get the offset


            }
        }
    }


}


__global__ void pushRelabel(deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp, densestCorePointer densestCore, deviceCliquesPointer finalCliqueData, ui *compCounter,ui *counter,ui *upperBound, ui *ranks,ui *offset,ui n, ui m, ui totalWarps, int totalComponents, ui k, ui lb) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(ui i = warpId; i < totalComponents; i += totalWarps){

        if(upperBound[warpId]>lb){
            ui alpha = (upperBound[warpId]+lb)/2;
            ui offsetLoc;
            ui start = conComp.componentOffset[i];
            ui end = conComp.componentOffset[i+1];
            ui total = end - start;
            ui cliqueStart = compCounter[i];
            ui cliqueEnd = compCounter[i+1];
            ui totalCliques = cliqueEnd-cliqueStart;

            ui rank = ranks[i];

            ui fStart = offset[rank];
            ui eStart = eOffset[rank];

            ui tFlow = totalCliques +  total + 2;

            ui s = tFlow - 2;
            ui t = tFlow - 1;

            for(ui i=laneId; i < tFlow; i+=warpsize){
                height[eStart + i ] = (i == s) ? total : 0;
                excess[eStart+ i] = 0;

            }
            __syncwarp();

            for(ui i=laneId; i < total; i+=warpsize){
                ui vertex = flowNetwork.toEdge[fStart + 2*totalCliques+i ];
                ui cap = flowNetwork.capacity[fStart + 2*totalCliques+i ];
                flowNetwork.flow[fStart + 2*totalCliques+i ] = cap;
                atomicAdd(&flowNetwork.excess[eStart + totalCliques+i ],cap);

            }

            __syncwarp();

            bool active = true;
            while (__any_sync(0xffffffff, active)){
                active = false;
                bool pushed = false;

                ui head = 2*totalCliques+ 3* total;

                // Push Vertex to sink
                for(ui i= head + laneId; i < total; i+=warpsize){
                    if (excess[eStart+i] > 0){

                        //Alawys Sink
                        ui v = flowNetwork.toEdge[offset + head + i ];
                        ui cap = flowNetwork.capacity[offset + head + i ];
                        ui flow = flowNetwork.flow[offset + head + i ];
                        ui resCap = cap-flow;
                        if(resCap>0 && flowNetwork.height[eStart+totalCliques+i]>flowNetwork.height[eStart+v]){
                            int push_amt = min(flowNetwork.excess[u], resCap);
                            atomicAdd(&flowNetwork.flow[offset + head + i], push_amt);
                            atomicSub(&flowNetwork.excess[eStart+ totalCliques + i], push_amt);
                            atomicAdd(&flowNetwork.excess[eStart+ totalCliques + v], push_amt);
                            pushed = true;
                            active = true;

                        }

                    }
                }
                
                //Check if we can push from cliques to verticies
                for(ui i=laneId; i < totalCliques; i+=warpsize){
                    if (excess[eStart+i] > 0){
                        
                        for(ui x = 0; x <k;x++){
                            ui v = flowNetwork.toEdge[offset + i* k +x ];
                            ui cap = flowNetwork.capacity[offset + i* k +x ];
                            ui flow = flowNetwork.flow[offset + i* k +x ];
                            ui resCap = cap-flow;

                            if(resCap>0 && flowNetwork.height[eStart+i]>flowNetwork.height[eStart+v]){
                                int push_amt = min(flowNetwork.excess[eStart+i], resCap);
                                atomicAdd(&flowNetwork.flow[offset + i* k +x ], push_amt);
                                atomicSub(&flowNetwork.excess[eStart+i], push_amt);
                                atomicAdd(&flowNetwork.excess[eStart+ v], push_amt);
                                pushed = true;
                                active = true;
                            }

                        }

                    }

                }

            
                for(ui i= laneId; i < tFlow; i+=warpsize){

                    int min_h = INT_MAX;
                    ui nStart = flowNetwork.offset[i];
                    for (ui j = flowNetwork.neighbors[nStart + i]; j < flowNetwork.neighbors[nStart + i+1]; ++j) {
                        int v = csr_col_ind[eid];
                        if (capacities[j] - flows[j] > 0)
                            min_h = min(min_h, height[v]);
                    }

                }
            
            }

            


        }
    }



}


__global__ void edmondsKarp(deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp, densestCorePointer densestCore, deviceCliquesPointer finalCliqueData, ui *compCounter,ui *counter,ui *upperBound, ui *ranks,ui *offset, int *augmentedPaths, ui *BFS, ui apSize, ui n, ui m, ui totalWarps, int totalComponents, ui k, ui lb){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    for(ui i = idx; i < totalComponents; i += TOTAL_THREAD){
        if(upperBound[i]>lb){
            ui alpha = (upperBound[i]+lb)/2;
            ui offsetLoc;
            ui start = conComp.componentOffset[i];
            ui end = conComp.componentOffset[i+1];
            ui total = end - start;

            ui bais = 1/(total*(total-1));

            ui cliqueStart = compCounter[i];
            ui cliqueEnd = compCounter[i+1];
            ui totalCliques = cliqueEnd-cliqueStart;

            ui offset = offset[rank[i]];

            ui apOffset = apSize*i;

            ui s = totalCliques + total;
            ui t = totalCliques + total + 1;

            ui curent = s;
            ui count = 1;
            ui tail = 0;
            BFS[apOffset] = s;
            ui maxFlow =INF;
            ui current;

            while( (u-l) > bais ){

                while(tail < count){
                    ui current = BFS[apOffset + tail];
                    tail++;
                    if(current == t){

                        while(current!=s){
                            ui nStart= flowNetwork.offset[offset + augmentedPaths[current]];
                            ui nEnd = flowNetwork.offset[offset + augmentedPaths[current]+1];
                            for(ui j = nStart; j < nEnd; j++ ){
                                if(j==current){
                                    ui availCap =  flowNetwork.capacity[nStart+j] -flowNetwork.flow[nStart+j];
                                    if(maxFlow > availCap){
                                        maxflow = availCap;
                                    }
                                    
                                    break;
                                }
                            }

                            current = augmentedPaths[current];


                        }

                    }

                    ui nStart= flowNetwork.offset[offset + current];
                    ui nEnd = flowNetwork.offset[offset + current+1];

                    for(ui j = nStart; j < nEnd; j++ ){
                        ui v = flowNetwork.toEdge[nStart+j];
                        ui availCap =  flowNetwork.capacity[nStart+j] -flowNetwork.flow[nStart+j];
                        if((augmentedPaths[apOffset+v] == -1) && availCap>0){
                            augmentedPaths[apOffset+v] = current;
                            BFS[apOffset+count] = v;
                            count++;
                        }
                    }

                }
                
            }


        }


    }
}


/*__global__ void createPaths(deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp, ui totalWarps, int totalComponents, ui k, ui alpha){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;
    ui offset = i*(2*totalCliques*k+4*total);
    for(ui i = warpId; i < totalComponents; i += totalWarps){
        



    }

}*/



/*__global__ void dynamicExact(deviceComponentPointers conComp,ui *newOffset,ui *newNeighbors, ui *components, ui n, ui m,ui totalComponents, ui totalWarps){

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;
    bool threadChanged = false;

    for(ui i = warpId; i < totalComponents; i += totalWarps){
        ui startComp = conComp.componentOffset[i];;
        ui endComp = conComp..componentOffset[i+1];
        ui total = endComp-startComp;



        //Create Flow network


        //get lower bound upper bound and bais

        //create spanning tree for forward edges

        //run algo

        // run algo again with backward edges 

        //atomic max for max density 



    }

}*/