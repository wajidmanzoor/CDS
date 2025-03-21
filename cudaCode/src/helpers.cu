#include "../utils/cuda_utils.cuh"
#include "../inc/helpers.cuh"

//TODO: remove unused variables 
//TODO: Make code a bit more structured and clean
//TODO: Try to find more ways to optimize

__device__ inline void writeToBuffer(ui *glBuffer, ui loc, ui v, ui bufferSize){
    assert(loc < bufferSize);
    glBuffer[loc] = v;
}

__device__ inline ui readFromBuffer(ui* glBuffer, ui loc,, ui bufferSize){
    assert(loc < bufferSize);
    return glBuffer[loc];
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
            //printf("warp %d lane %i vertex %d neighoff %d canoff %d neigh %d label %d \n",i,j,vertex, neighOffset,candidateOffset,neigh,label[warpId*n + neigh]);
            if(label[warpId*n + neigh] == k) {
                label[warpId*n + neigh] = k - 1;
                //printf("WarpId %d laneId %d candidate %d num neigh %d neigh %d label %d \n", warpId,laneId,vertex,D.degree[vertex],neigh,label[warpId*n + neigh]);
                ui loc = atomicAdd(&counter[threadIdx.x / warpSize], 1);
                levelData.candidatesPartition[candidateOffset + loc] = neigh;
            }
        }
        __syncwarp();
        if(laneId == 0 && counter[threadIdx.x / warpSize] > 0) {
            //printf("warp %d counter %d cliquePart %d leveldataC %d vertex %d \n",i,counter[threadIdx.x / warpSize], cliquePartition, levelData.count[warpId + 1], vertex);
            levelData.partialCliquesPartition[cliquePartition + levelData.count[warpId + 1] * (k-1) + level] = vertex;
            levelData.count[warpId + 1] += 1;
            levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]] =
                levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1] - 1] + counter[threadIdx.x / warpSize];
        }
         __syncwarp();

        int start = candidateOffset;
        //int end = candidateOffset + counter[threadIdx.x / warpSize];
        //printf("warp %d lane Id %d start %d end %d total %d \n", warpId, laneId, start, end, counter[threadIdx.x / warpSize]);
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

                //printf("WarpId %d laneId %d candidate %d num neigh %d   mask part %d count %d offset loc %d offset %d max bit %d  bitIndex %d loc %d  bitmask %d \n", warpId,laneId,candidate,degree,maskPartition,levelData.count[warpId + 1],offsetPartition + levelData.count[warpId + 1] -1,levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1] -1 ],maxBitMask,bitmaskIndex,maskPartition + (levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]-1]+j) * maxBitMask + bitmaskIndex,bitmask);


                levelData.validNeighMaskPartition[maskPartition + (levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]-1]+j) * maxBitMask + bitmaskIndex] = bitmask;
            }
        }

        __syncwarp();

        for(int i = laneId; i<n;i+=32){
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
            
            int totalMasks = (D.degree[candidate]+31)/32;
            for(int j =0; j < totalMasks; j++){
                levelData.validNeighMask[(writeOffset+i)*maxBitMask + j ] =
                levelData.validNeighMaskPartition[maskPartition + (levelData.offsetPartition[offsetPartition + iter] + i)*maxBitMask + j];
            }

        }

        if(laneId< level+1 ){

                levelData.partialCliques[levelData.count[warpId]*(k-1)+ iter*(k-1) + laneId] = levelData.partialCliquesPartition[cliquePartition + iter * (k-1) + laneId];
                //printf("warp Id %d lane id %d count %d iter % d pclique %d level %d write loc %d \n",warpId,i,levelData.count[warpId],iter,levelData.partialCliquesPartition[cliquePartition + iter * (k-1) + i],i,levelData.count[warpId]+ iter*(k-1) + i);
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
        //printf("warp %d start %d total %d\n", warpId, start, totalCandidates);

        for(int iter = 0; iter <totalCandidates; iter ++){
            int candidate = levelData.candidates[start + iter];
            if(laneId==0){
                counter[threadIdx.x/warpSize] = 0;
            }
            __syncwarp();
            //printf("---warp %d start %d iter %d cand %d \n", warpId, start, iter,candidate);

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
                //printf("warp %d lane id %d start %d iter %d cand %d neigh %d label %d loc %d \n", warpId, j, start, iter,candidate,neigh,label[warpId*n + neigh],writeOffset + loc);

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
    
                    //printf("WarpId %d laneId %d candidate %d num neigh %d   mask part %d count %d offset loc %d offset %d max bit %d  bitIndex %d loc %d  bitmask %d \n", warpId,laneId,candidate,degree,maskPartition,levelData.count[warpId + 1],offsetPartition + levelData.count[warpId + 1] -1,levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1] -1 ],maxBitMask,bitmaskIndex,maskPartition + (levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]-1]+j) * maxBitMask + bitmaskIndex,bitmask);
    
    
                    levelData.validNeighMaskPartition[maskPartition + (levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]-1]+j) * maxBitMask + bitmaskIndex] = bitmask;
                }
            }
    
            __syncwarp();
    
            for(int i = laneId; i<n;i+=32){
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
        //printf("warp %d start %d total %d\n", warpId, start, totalCandidates);

        for(int iter = 0; iter <totalCandidates; iter ++){
            int candidate = levelData.candidates[start + iter];
            if(laneId==0){
                counter[warpId]=0;
            }
            int degree = D.degree[candidate];
            int neighOffset = D.offset[candidate];
            //printf("---warp %d start %d iter %d cand %d \n", warpId, start, iter,candidate);

            
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

        printf("idx %d \n", idx);
        for(int j=0;j<k;j++){
            
            elements[j] = cliqueData.trie[j*t+i];
            degree[j] = G.cliqueDegree[elements[j]];
            //printf("idx %d j %d loc element %d element %d degree %d \n",idx,j,j*t+i,elements[j],degree[j]);
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
            //printf("idx %d j %d loc element %d element %d \n",idx,j,j*t+i,elements[j]);
    
        }

    }

}


__global__ void selectNodes(deviceGraphPointers G, ui *bufTails,ui *glBuffers, ui glBufferSize, ui n, ui level){
    __shared__ ui *glBuffer; 
    __shared__ ui bufTail; 
    
    if(THID == 0){
        bufTail = 0;
        glBuffer = glBuffers + blockIdx.x*glBufferSize;
    }
    __syncthreads();

    ui idx = blockIdx.x * blockDim.x + threadIdx.x; 
    for(ui base = 0; base < n; base += BLK_DIM){
        
        ui v = base + idx; 

        if(v >= n) continue;

        if(G.cliqueCore[v] == level){
            ui loc = atomicAdd(&bufTail, 1);
            writeToBuffer(glBuffer, loc, v,glBufferSize);
        }
    }

    __syncthreads();

    if(THID == 0) 
    {
        bufTails [blockIdx.x] = bufTail;
    }

}

__global__ void processNodesByWarp(deviceGraphPointers G,deviceCliquesPointer cliqueData, ui *bufTails,ui *glBuffers, ui *globalCount, ui glBufferSize, ui n, ui level, ui k){
    __shared__ ui bufTail;
    __shared__ ui *glBuffer;
    __shared__ ui base;
    ui warpId = THID / 32;
    ui laneId = THID % 32;
    ui regTail;
    ui i;
    if(THID==0){
    bufTail = bufTails[blockIdx.x];
    base = 0;
    glBuffer = glBuffers + blockIdx.x*GLBUFFER_SIZE; 
    assert(glBuffer!=NULL);
    }

    while(true){
    __syncthreads(); 
    if(base == bufTail) break; // all the threads will evaluate to true at same iteration
    i = base + warp_id;
    regTail = bufTail;
    __syncthreads();

    if(i >= regTail) continue; // this warp won't have to do anything            

    if(THID == 0){
    base += WARPS_EACH_BLK;
    if(regTail < base )
    base = regTail;
    }
    //bufTail is incremented in the code below:
    ui v = readFromBuffer(glBuffer, i,bufferSize);

 
    __syncwarp();
    for(ui i =laneId; i<t; i+=warpSize){

        if( (v = cliqueData.trie[i]) && (cliqueData.status[i]==1)){
            for(ui j =1;j<k;j++){
                ui u = cliqueData.trie[j*t+i];
                int a = atomicSub(&G.cliqueCore[u]);
                if(a == level+1){
                    ui loc = atomicAdd(&bufTail, 1);
            
                    writeToBuffer(glBuffer, loc, u);
                }
                if(a <= level){ 
                    atomicAdd(&G.cliqueCore[u], 1);
                }
            }
            cliqueData.status[i]==0;
            
            
        }
    }
    

    if(THID == 0 && bufTail>0){
    atomicAdd(globalCount, bufTail); // atomic since contention among blocks
    }
}
}

__global__ void processNodesByBlock(deviceGraphPointers G,deviceCliquesPointer cliqueData, ui *bufTails,ui *glBuffers, ui *globalCount, ui glBufferSize, ui n, ui level, ui k){
    __shared__ ui bufTail;
    __shared__ ui *glBuffer;
    __shared__ ui base;
    ui warpId = THID / 32;
    ui laneId = THID % 32;
    ui regTail;
    ui i;
    if(THID==0){
    bufTail = bufTails[blockIdx.x];
    base = 0;
    glBuffer = glBuffers + blockIdx.x*GLBUFFER_SIZE; 
    assert(glBuffer!=NULL);
    }

    while(true){
    __syncthreads(); 
    if(base == bufTail) break; // all the threads will evaluate to true at same iteration
    i = base + blockIdx.x;
    regTail = bufTail;
    __syncthreads();

    if(i >= regTail) continue; // this warp won't have to do anything            

    if(THID == 0){
    base += 1;
    if(regTail < base )
    base = regTail;
    }
    //bufTail is incremented in the code below:
    ui v = readFromBuffer(glBuffer, i,bufferSize);

 
    __syncwarp();
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    for(ui i = idx; i<t; i+= BLK_DIM){
        if( (v = cliqueData.trie[i]) && (cliqueData.status[i]==1)){
            for(ui j =1;j<k;j++){
                ui u = cliqueData.trie[j*t+i];
                int a = atomicSub(&G.cliqueCore[u]);
                if(a == level+1){
                    ui loc = atomicAdd(&bufTail, 1);
            
                    writeToBuffer(glBuffer, loc, u);
                }
                if(a <= level){ 
                    atomicAdd(&G.cliqueCore[u], 1);
                }
            }
            cliqueData.status[i]==0;
            
            
        }

    }
    

    if(THID == 0 && bufTail>0){
    atomicAdd(globalCount, bufTail); // atomic since contention among blocks
    }
}
}

