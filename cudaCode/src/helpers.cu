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


    // Each Warp processes on vertex
    for(int i = warpId; i < n; i += totalWarps) {

        ui vertex = i;
        ui neighOffset = D.offset[vertex];
        if(laneId == 0) {
            counter[threadIdx.x / warpSize] = 0;
        }

        __syncwarp();

        int candidateOffset = candidatePartition + levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]];

        // Iter through neighbors of vertex and add to candidate in label = k.
        for(int j = laneId; j < D.degree[vertex]; j += warpSize) {
            ui neigh = D.neighbors[neighOffset + j];
            //printf("warp %d lane %i vertex %d neighoff %d canoff %d neigh %d label %d \n",i,j,vertex, neighOffset,candidateOffset,neigh,label[warpId*n + neigh]);

            ui old_val = atomicCAS(&label[warpId*n + neigh], k, k-1);
            if(old_val==k) {
                //printf("WarpId %d laneId %d candidate %d num neigh %d neigh %d label %d \n", warpId,laneId,vertex,D.degree[vertex],neigh,label[warpId*n + neigh]);
                ui loc = atomicAdd(&counter[threadIdx.x / warpSize], 1);
                //assert((candidateOffset + loc) > ((warpId+1) * cpSize));
                levelData.candidatesPartition[candidateOffset + loc] = neigh;
            }
        }
        __syncwarp();

        // Offset and PC added
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

        for(int x = laneId; x<n;x+=32){
          label[warpId*n + x] = k;
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
            levelData.candidates[writeOffset+ i] = candidate;

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

    // Iter through total pc.
    for(int i =warpId; i < totalTasks ; i+= totalWarps ){

        int start = levelData.offset[i];
        int totalCandidates = levelData.offset[i+1]- start;
        //printf("total warp %d start %d total %d\n", warpId, start, totalCandidates);

        // Process one canddate at a time
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
                int neighBitMask = levelData.validNeighMask[(start+iter)*maxBitMask + iterBitMask];
                ui neigh = D.neighbors[neighOffset + j];
                if(neighBitMask & (1 << bitPos )){
                     ui old_val = atomicCAS(&label[warpId*n + neigh], iterK, iterK-1);
                     if(old_val== iterK){
                        ui loc = atomicAdd( &counter[threadIdx.x/warpSize], 1);
                        levelData.candidatesPartition[writeOffset + loc] = neigh;

                     }
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


                        if (label[warpId*n + D.neighbors[neighOffset + x]] == (iterK - 1)) {
                            bitmask |= (1 << (x - startNeighbor)); // Set the bit for valid neighbors
                            //printf("Warpid %d lane %d iter %d can %d x %d negh %d label %d bit mask %d write loc %d ind %d \n",warpId,laneId,iter,cand,x,D.neighbors[neighOffset + x],label[warpId*n + D.neighbors[neighOffset + x]],bitmask,(levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]-1]+j),bitmaskIndex);


                        }
                    }

                    levelData.validNeighMaskPartition[maskPartition + (levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]-1]+j) * maxBitMask + bitmaskIndex] = bitmask;
                }
            }

            __syncwarp();

            for(int x = laneId; x<n;x+=32){
              label[warpId*n + x] = iterK;
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
                counter[threadIdx.x/warpSize]=0;
            }
            __syncwarp();
            int degree = D.degree[candidate];
            int neighOffset = D.offset[candidate];
            //printf("---warp %d start %d iter %d cand %d \n", warpId, start, iter,candidate);


            for(int j = laneId; j< degree; j+= warpSize ){
                int iterBitMask = j/warpSize;
                int bitPos = j%32;
                int neighBitMask = levelData.validNeighMask[(start+iter)*maxBitMask + iterBitMask];
                if(neighBitMask & (1 << bitPos )){

                    ui neigh = D.neighbors[neighOffset + j];

                    ui loc = atomicAdd(globalCounter,1);
                    for(int ind =0; ind < k-2; ind++){
                        cliqueData.trie[trieSize * ind + loc] = levelData.partialCliques[(i)*(k-1) + ind];

                    }
                    atomicAdd(&counter[threadIdx.x/warpSize],1);
                    cliqueData.trie[trieSize * (k-2) + loc]  = candidate;
                    cliqueData.trie[trieSize * (k-1) + loc] = neigh;
                    cliqueData.status[loc]= -1;
                    atomicAdd(&G.cliqueDegree[neigh],1);
                    atomicAdd(&G.cliqueDegree[candidate],1);

                }


            }
            __syncwarp();

            for(int j = laneId; j< k-2 ; j+= warpSize ){
                int pClique = levelData.partialCliques[i*(k-1) + j];
                atomicAdd(&G.cliqueDegree[pClique],counter[threadIdx.x/warpSize]);
            }
             __syncwarp();

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

        //printf("idx %d \n", idx);
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

__global__ void selectNodes(deviceGraphPointers G, ui *bufTails, ui *glBuffers, ui glBufferSize, ui n, ui level) {
    __shared__ ui *glBuffer;
    __shared__ ui bufTail;

    if(threadIdx.x == 0) {
        bufTail = 0;
        glBuffer = glBuffers + blockIdx.x * glBufferSize;
    }
    __syncthreads();
    ui total = (n+BLK_NUMS)/BLK_NUMS;
    ui start = blockIdx.x * total;
    ui end = min(start + total, n);

    ui thid = threadIdx.x;
    /*if(thid ==0){
        printf("blockid %d start %d end %d \n",blockIdx.x,start,end);
    }*/

    for(ui  i = start + thid; i<end ; i+= BLK_DIM) {
        ui v = i;
        if(G.cliqueCore[v] == level) {
            ui loc = atomicAdd(&bufTail, 1);
            assert(loc < glBufferSize && "glBuffer overflow in selectNodes");
            glBuffer[loc] = v;
        }
    }

    __syncthreads();

    if(threadIdx.x == 0) {
        assert(bufTail <= glBufferSize && "bufTail exceeds glBufferSize in selectNodes");
        bufTails[blockIdx.x] = bufTail;
    }
}

__global__ void processNodesByWarp(deviceGraphPointers G,deviceCliquesPointer cliqueData, ui *bufTails,ui *glBuffers, ui *globalCount, ui glBufferSize, ui n, ui level, ui k, ui totalCliques){
    /*removes the verties to get their core value.
      Warp processes the verticies in its virtual partition parallely*/
    __shared__ ui bufTail;
    __shared__ ui *glBuffer;
    __shared__ ui base;
    ui warpId = threadIdx.x / 32;
    ui laneId = threadIdx.x % 32;
    ui regTail;
    ui i;
    if(threadIdx.x==0){
        // index of last vertex in vitual partition of glguffer for current warp
        bufTail = bufTails[blockIdx.x];

        // stores index of current processed vertex
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
            if(regTail < base );
            base = regTail;
        }
        //__syncthreads();


        //vertex to be removed
        ui v = glBuffer[i];


        __syncwarp();

        // warp iters through the clique data.
        for(ui j =laneId; j< totalCliques; j+=warpSize){

            // if valid clique and not removed yet
            if(cliqueData.status[j] == -1){

                // flag to check if vertex found in the clique
                bool found = false;
                // stores the index at which the vertex was found in clique (0,k-1)
                ui w =0;

                // iter through verticies of clique sequentially.
                while(w<k){

                    if(cliqueData.trie[w*totalCliques+j] == v){

                    found = true;
                    //break;
                    }
                    w++;
                }

                if(found){
                    cliqueData.status[j] = level;
                    // iter throught the clique verticies
                    for(ui x =0 ;x<k;x++){

                        // continue it clique vertex is same as vertex
                        // clique vertex
                        ui u = cliqueData.trie[x*totalCliques+j];
                        if(u==v) continue;

                        // decreament its core value by 1
                        int a = atomicSub(&G.cliqueCore[u],1);

                        // if core value is less than level, update to level.
                        if (a-1 < level) {
                            atomicExch(&G.cliqueCore[u], level);
                        }
                        // if core value is level, add to glbuffer so can be removed in this level.
                        if(a == level+1){
                            ui loc = atomicAdd(&bufTail, 1);
                            glBuffer[loc] = u;

                        }
                    }

                    // set status of the clique to current core level.



                }
            }
        }
    }

    __syncthreads();
    if(threadIdx.x==0){
        atomicAdd(globalCount, bufTail); // atomic since contention among blocks
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
        //printf("blockId %d idx %d base %d i %d reg tail %d vertex %d cc %u status %d \n", blockIdx.x,idx,base,i,regTail,v,cliqueData.trie[j],cliqueData.status[j]);

        if( (v == cliqueData.trie[j]) && (cliqueData.status[j] == -1 )){
            for(ui x =1;x<k;x++){
                ui u = cliqueData.trie[x*t+j];
                int a = atomicSub(&G.cliqueCore[u], 1);
                if(a == level+1){
                    ui loc = atomicAdd(&bufTail, 1);
                    glBuffer[loc] = u;

                }
                if(a <= level){
                    atomicAdd(&G.cliqueCore[u], 1);
                }
                if(G.cliqueCore[u]<0){
                  G.cliqueCore[u] = 0;
                }
            }
            cliqueData.status[j] = level;


        }

    }


    __syncthreads();

    if(threadIdx.x == 0 && bufTail>0){
    atomicAdd(globalCount, 1); // atomic since contention among blocks
    }
}
}

__global__ void generateDensestCore(deviceGraphPointers G, densestCorePointer densestCore,ui *globalCount, ui n, ui core, ui totalWarps){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(ui i = warpId; i < n; i += totalWarps){
        if(G.cliqueCore[i]>= core){

            //printf("warpid %d lane %d core %d density %d i %d \n",warpId,laneId,G.cliqueCore[i],density,i);
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
                if(G.cliqueCore[neigh] >= core) {
                    count++;
                }
            }
            for (int offset = warpSize / 2; offset > 0; offset /= 2) {
                count += __shfl_down_sync(0xFFFFFFFF, count, offset);
            }
            __syncwarp();

            if(laneId == 0) {
                densestCore.offset[loc+1] = count;
            }
            __syncwarp();


        }
    }
}

__global__ void generateNeighborDensestCore(deviceGraphPointers G, densestCorePointer densestCore, ui core, ui totalWarps) {

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

            if(G.cliqueCore[neigh] >= core) {
                int loc = atomicAdd(&counter[threadIdx.x / warpSize], 1);

                densestCore.neighbors[loc] = densestCore.reverseMap[neigh];

            }
        }
      __syncwarp();
    }
}

__global__ void pruneEdges(densestCorePointer densestCore, deviceCliquesPointer cliqueData, ui *pruneStatus,ui totalCliques, ui k, ui level ){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(ui i = warpId; i<totalCliques; i+=TOTAL_WARPS){

        if(cliqueData.status[i] >= level){

            for(ui iter =0; iter< k ; iter ++){
                // v should be mapped
                ui u_ = ((iter)%k)*totalCliques+i;
                ui u  =  densestCore.reverseMap[cliqueData.trie[u_]];
                for(ui j = 0; j < k; j++){
                    ui v_ = ((j)%k)*totalCliques+i;


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

__global__ void getConnectedComponentStatus(deviceComponentPointers conComp,deviceCliquesPointer cliqueData, densestCorePointer densestCore, ui *compCounter, ui totalCliques, ui k,ui maxCore, ui totalThreads){

    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    for(ui i =idx; i<totalCliques;i +=totalThreads){
        if(cliqueData.status[i]>=maxCore){
          int comp = INT_MAX;

          for(ui x=0;x<k;x++){
            ui vertex =densestCore.reverseMap[cliqueData.trie[x*totalCliques+ i]];
            comp = min(comp,conComp.components[conComp.reverseMapping[vertex]]);

            cliqueData.trie[x*totalCliques + i] = vertex;


          }
          cliqueData.status[i] = comp;
          atomicAdd(&compCounter[comp+1],1);


      }else{
        cliqueData.status[i] = -1;
      }
    }

}

__global__ void rearrangeCliqueData(deviceComponentPointers conComp,deviceCliquesPointer cliqueData, deviceCliquesPointer finalCliqueData,densestCorePointer densestCore, ui *compCounter,ui *counter,ui totaLCliques,ui k, ui newTotaLCliques ,ui totalThreads){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    for(ui i =idx; i<totaLCliques;i +=totalThreads){

        int comp =  cliqueData.status[i];

        if(comp>-1){
          ui loc;
          for(ui j=0;j<k;j++){
              ui vertex = cliqueData.trie[j*totaLCliques + i];
              ui offset = compCounter[comp];
              if(j==0){
                  loc = atomicAdd(&counter[comp],1);
                  finalCliqueData.status[offset + loc ] = comp;


              }
              finalCliqueData.trie[offset + j*newTotaLCliques + loc ] = vertex;
          }
        }

    }

}

__global__ void getLbUbandSize(deviceComponentPointers conComp, ui *compCounter, double *lowerBound, double *upperBound, ui *ccOffset,  ui *neighborSize, ui totalComponenets, ui k, double maxDensity){
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if(idx==0){
        ccOffset[idx]=0;
        neighborSize[idx]=0;

    }

    for(ui i = idx; i<totalComponenets; i+=TOTAL_THREAD ){
        ui totalCliques = compCounter[i+1] - compCounter[i];
        ui totalSize = conComp.componentOffset[i+1] -  conComp.componentOffset[i];
        double lb = (double) (totalCliques)/totalSize;
        lowerBound[i]  = lb;

        double dem = pow(fact(k),1.0/k);
        double num = pow(totalCliques, (k-1.0)/k);
        double ub = min(maxDensity, num/dem);

        upperBound[i] = ub;

        if(ub>lb){
            ccOffset[i+1] = totalCliques + totalSize + 2 +1 ;
            neighborSize[i+1] = 2*(2*totalCliques*k + 2*totalSize)+1;


        }
        else{
            ccOffset[i+1] = 0;
            neighborSize[i+1] = 0;

        }
    }


}


__global__ void createFlowNetworkOffset(deviceGraphPointers G, deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp, deviceCliquesPointer finalCliqueData, ui *compCounter,double *upperBound , ui totalWarps, ui totalComponents, ui k, double lb, ui t){



    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(ui i = warpId; i < totalComponents; i += totalWarps){

        if(upperBound[i]>lb){

            ui start = conComp.componentOffset[i];
            ui end = conComp.componentOffset[i+1];
            ui total = end - start;
            ui startClique = compCounter[i];
            ui totalCliques = compCounter[i+1]-compCounter[i];

            ui vertexOffset = flowNetwork.offset[i];
            //ui neighborOffset = flowNetwork.neighborOffset1[i];

            //printf("warpid %d laneId %d start %d end %d total %d neighborOffset %d \n",warpId,laneId,start,end,total,neighborOffset);

            for (ui j = laneId; j < total; j += warpSize){
                ui vertex = conComp.mapping[start+j];

                ui cliqueDegree = 0;
                for(ui x =0; x < totalCliques; x ++){
                  ui u;
                  for(ui k_ = 0; k_<k;k_++){
                    u = finalCliqueData.trie[t*k_ + startClique+x];
                    if (u==vertex){
                      cliqueDegree++;
                      }
                  }


                }
                flowNetwork.neighborOffset2[vertexOffset+j + 1] = 2*(cliqueDegree + 1);
            }
            for (ui j = laneId; j < totalCliques; j += warpSize){
                flowNetwork.neighborOffset2[vertexOffset+total+j+1] = 2*k;
            }
            if(laneId==0){
                flowNetwork.neighborOffset2[vertexOffset+total+totalCliques+1] = total;
                flowNetwork.neighborOffset2[vertexOffset+total+totalCliques+2] = total;
                flowNetwork.neighborOffset2[0]=0;
            }


            }
        }
    }

__global__ void createFlowNetwork(deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp, deviceCliquesPointer finalCliqueData, ui *compCounter,double *upperBound , ui totalWarps, ui totalComponents, ui k, double lb, ui t){



    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(ui i = warpId; i < totalComponents; i += totalWarps){

        if(upperBound[i]>lb){

            ui start = conComp.componentOffset[i];
            ui end = conComp.componentOffset[i+1];
            ui total = end - start;
            ui startClique = compCounter[i];
            ui totalCliques = compCounter[i+1]-compCounter[i];

            ui vertexOffset = flowNetwork.offset[i];

            double alpha = (upperBound[i]+lb)/2;

            //printf("warpid %d laneId %d start %d end %d total %d neighborOffset %d \n",warpId,laneId,start,end,total,neighborOffset);

            for (ui j = laneId; j < total; j += warpSize){
                ui neighborOffset = flowNetwork.neighborOffset2[vertexOffset+j];


                // Vertex to sink
                flowNetwork.Edges[neighborOffset] = total+totalCliques+1;
                flowNetwork.capacity[neighborOffset] = alpha * k;

                //vertex to source (backward)
                flowNetwork.Edges[neighborOffset+1] = total+totalCliques;
                flowNetwork.capacity[neighborOffset+1] = 0.0;


                ui vertex = conComp.mapping[start+j];

                ui cliqueDegree = flowNetwork.neighborOffset2[vertexOffset+j + 1] - flowNetwork.neighborOffset2[vertexOffset+j];


                ui temp= 2;

                for(ui x =0; x < totalCliques; x ++){
                    ui u;
                    for(ui k_ = 0; k_<k;k_++){
                      u = finalCliqueData.trie[t*k_ + startClique+x];
                      if (u==vertex){
                        //vertex to clique
                        flowNetwork.Edges[neighborOffset+ temp] = total+x;
                        flowNetwork.capacity[neighborOffset+ temp] = 1.0;

                        // vertex to clique backward
                        flowNetwork.Edges[neighborOffset+ temp+ 1] = total+x;
                        flowNetwork.capacity[neighborOffset+ temp+1] = 0.0;
                        temp+=2;
                        }
                    }
                    if(temp==(2*cliqueDegree)){
                        break;
                    }


                  }

            }
            for (ui j = laneId; j < totalCliques; j += warpSize){
                ui neighborOffset = flowNetwork.neighborOffset2[vertexOffset+total+j];
                ui u;
                for(ui k_ = 0; k_<k;k_++){
                    u = finalCliqueData.trie[t*k_ + startClique+j];

                    //Clique to vertex
                    flowNetwork.Edges[neighborOffset+ 2*k_ ] = conComp.reverseMapping[u] -start;
                    flowNetwork.capacity[neighborOffset+ 2*k_] = DINF;

                    //Clique to vertex backward
                    flowNetwork.Edges[neighborOffset+ 2*k_ +1 ] = conComp.reverseMapping[u] -start;
                    flowNetwork.capacity[neighborOffset+ 2*k_ + 1] = 0;

                }


            }
            ui neighborOffset_source = flowNetwork.neighborOffset2[vertexOffset+total+totalCliques];
            ui neighborOffset_sink = flowNetwork.neighborOffset2[vertexOffset+total+totalCliques +1 ];



            for (ui j = laneId; j < total; j += warpSize){

                ui cliqueDegree = (flowNetwork.neighborOffset2[vertexOffset+j + 1] - flowNetwork.neighborOffset2[vertexOffset+j]-2)/2;

                //source to vertex
                flowNetwork.Edges[neighborOffset_source+ j] = j;
                flowNetwork.capacity[neighborOffset_source+ j] = (double) cliqueDegree;


                //sink to vertex backward
                 flowNetwork.Edges[neighborOffset_sink+ j] = j;
                flowNetwork.capacity[neighborOffset_sink+ j] = 0.0;


            }



            }

    }
 }

__global__ void pushRelabel(deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp,  deviceCliquesPointer finalCliqueData, ui * compCounter, double * upperBound, double * lowerBound, ui * activeNodes, ui * componenetsLeft, ui * checkResult, ui totalWarps, int totalComponents, ui k, ui t, ui partitionSize) {
    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui * counter = (ui * )(sharedMemory + sizeOffset);

    // += WARPS_EACH_BLK * sizeof(ui);
    //sizeOffset = (sizeOffset + alignof(double) - 1) & ~(alignof(double) - 1);
    //ui *densities = (ui *)(sharedMemory + sizeOffset);

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for (ui i = warpId; i < totalComponents; i += totalWarps) {
        ui start = conComp.componentOffset[i];
        ui end = conComp.componentOffset[i + 1];
        ui total = end - start;

        ui cliqueStart = compCounter[i];
        ui cliqueEnd = compCounter[i + 1];
        ui totalCliques = cliqueEnd - cliqueStart;

        ui fStart = flowNetwork.offset[i];

        ui tFlow = totalCliques + total + 2;

        double bais = 1.0 / (tFlow * (tFlow - 1));
        if (bais < 0.0001) {
        bais = 0.0001;
        }

        if ((upperBound[i] - lowerBound[i]) > bais) {

        ui s = tFlow - 2;
        ui t = tFlow - 1;

        //Set Height to 0 expect for s to total
        for (ui j = laneId; j < tFlow; j += warpSize) {
            flowNetwork.height[fStart + j] = (j == s) ? tFlow : 0;
            flowNetwork.excess[fStart + j] = 0;

        }
        __syncwarp();

        // Send Intial flow for s to all v
        for (ui j = laneId; j < total; j += warpSize) {
            ui nStart = flowNetwork.neighborOffset2[fStart + totalCliques + total];
            ui neigh = flowNetwork.Edges[nStart + j];
            double cap = flowNetwork.capacity[nStart + j];

            //Forward Flow s to vertex
            flowNetwork.flow[nStart + j] = cap;
            atomicAdd( & flowNetwork.excess[fStart + neigh], cap);

            //BackwardFlow vertex to s
            flowNetwork.flow[flowNetwork.neighborOffset2[fStart + neigh] + 1] = -cap;

        }

        __syncwarp();

        int maxIterations = 1000;

        //printf("warpid %d laneid %d\n ", warpId,laneId);

        //Push or Relabel until converge
        for (int iter = 0; iter < maxIterations; iter++) {
            if (laneId == 0) {
                counter[threadIdx.x / warpSize] = 0;
                //printf("reset iter %d warp %d c %d\n", iter, warpId, counter[threadIdx.x / warpSize]);
            }

            __syncwarp();


            //check for nodes with excess
            for (ui j = laneId; j < tFlow; j += warpSize) {
                if ((j != s) && (j != t) && (flowNetwork.excess[fStart + j] > 1e-10)) {
                    int pos = atomicAdd( & counter[threadIdx.x / warpSize], 1);
                    activeNodes[i * partitionSize + pos] = j;
                }
            }

            __syncwarp(); // Ensure all threads have the latest counter value

            bool should_break = (counter[threadIdx.x / warpSize] == 0);
            //printf("iter %d warpId %d laneid %d c %d  should_break %d \n", iter, warpId, laneId,counter[threadIdx.x / warpSize],should_break);

            unsigned mask = __ballot_sync(0xFFFFFFFF, should_break);
            if (mask == 0xFFFFFFFF) {  // All threads agree we should break

                //printf("break iter %d warpId %d lane %d  c %d \n", iter, warpId,laneId, counter[threadIdx.x / warpSize]);

                break;
            }

            bool pushed = false;

            for (ui j = laneId; j < counter[threadIdx.x / warpSize]; j += warpSize) {
                ui vertex = activeNodes[i * partitionSize + j];

                ui nStart = flowNetwork.neighborOffset2[fStart + vertex];
                ui nEnd = flowNetwork.neighborOffset2[fStart + vertex + 1];

                //Check neighbors to send acess to.
                for (ui x = nStart; x < nEnd; x++) {
                    //
                    /*if (flowNetwork.excess[fStart + vertex] == 0) {
                    break;
                    }*/
                    ui neigh = flowNetwork.Edges[x];
                    double residual = flowNetwork.capacity[x] - flowNetwork.flow[x];

                    // __syncwarp();

                    // If neighbor has capacity
                    if ((flowNetwork.height[fStart + vertex] == flowNetwork.height[fStart + neigh] + 1) && residual > 0) {

                        //printf("warp Id %d laneid %d vertex %d nstart %d end %d neigh %d resi %f cap %f flow %f loc %d hv %d hu %d ev %f en %f \n", i, j, vertex, nStart, nEnd, neigh, residual, flowNetwork.capacity[x], flowNetwork.flow[x], x, flowNetwork.height[fStart + vertex], flowNetwork.height[fStart + neigh], flowNetwork.excess[fStart + vertex], flowNetwork.excess[fStart + neigh]);

                        double delta = min(flowNetwork.excess[fStart + vertex], residual);
                        if (delta > 1e-10) {

                            //forward flow vertex to neigh
                            atomicAdd( & flowNetwork.flow[x], delta);

                            //Backward Flow neigh to vertex
                            ui stemp = flowNetwork.neighborOffset2[fStart + neigh];
                            ui etemp = flowNetwork.neighborOffset2[fStart + neigh + 1];
                            for (ui ind = stemp; ind < etemp; ind++) {
                                if (flowNetwork.Edges[ind] == vertex) {
                                    atomicAdd( & flowNetwork.flow[ind], -delta);
                                    //break;
                                }
                            }

                            //Decrease vertex excess
                            atomicAdd( & flowNetwork.excess[fStart + vertex], -delta);

                            //Increase neighbor excess
                            atomicAdd( & flowNetwork.excess[fStart + neigh], delta);
                            pushed = true;

                        }
                        // printf("after -- warp Id %d laneid %d vertex %d nstart %d end %d neigh %d resi %f cap %f flow %f loc %d hv %d hu %d ev %f en %f \n", i, j, vertex, nStart, nEnd, neigh, residual, flowNetwork.capacity[x], flowNetwork.flow[x], x, flowNetwork.height[fStart + vertex], flowNetwork.height[fStart + neigh], flowNetwork.excess[fStart + vertex], flowNetwork.excess[fStart + neigh]);

                    }

                }

                //__syncwarp();

                // Relabel
                if (!pushed && flowNetwork.excess[fStart + vertex] > 0) {
                    ui minHeight = UINT_MAX;
                    //printf("Relabel  WarpId %d laneid %d vertex %d excess %f height %d nstart %d end %d \n", i, laneId, vertex, flowNetwork.excess[fStart + vertex], flowNetwork.height[fStart + vertex], nStart, nEnd);
                    for (ui x = nStart; x < nEnd; x++) {
                    ui neigh = flowNetwork.Edges[x];
                    double residual = flowNetwork.capacity[x] - flowNetwork.flow[x];
                    //printf("warp %d laneid %d j %d v %d res %f negh %d height %d nstart %d \n", i, laneId, j, vertex, residual, neigh, flowNetwork.height[fStart + neigh], nStart);

                    if (residual > 1e-10) {
                        minHeight = (flowNetwork.height[fStart + neigh] < minHeight) ? flowNetwork.height[fStart + neigh] : minHeight;
                    }

                    }
                    if (minHeight != UINT_MAX) {
                    flowNetwork.height[fStart + vertex] = minHeight + 1;

                    }

                    //printf("after Relabel  WarpId %d laneid %d vertex %d excess %f height %d minHeight %d \n", i, laneId, vertex, flowNetwork.excess[fStart + vertex], flowNetwork.height[fStart + vertex], minHeight);

                }

            }

            __syncwarp();




        }

        __syncwarp();


        //printf("second warpid %d laneid %d\n ", warpId,laneId);

        if (laneId == 0) {
            //getRes[threadIdx.x / warpSize] = 0;
            double alpha = (upperBound[i] + lowerBound[i]) / 2;
            double temp = totalCliques * k;

            //printf("Warp %d lane %d excess %f temp %f \n",warpId,laneId,flowNetwork.excess[fStart + t],temp);
             double epsilon = 1e-3;
            if (fabs(flowNetwork.excess[fStart + t] - temp) < epsilon) {
              upperBound[i] = alpha;
              //printf("up here");
            } else {
              lowerBound[i] = alpha;
              checkResult[i] = 1;
              //printf("lb here");

            }
            if ((upperBound[i] - lowerBound[i]) > bais) {
            atomicAdd(componenetsLeft, 1);

            }

        }

        __syncwarp();

        //Update the network
        if ((upperBound[i] - lowerBound[i]) > bais){
                double alpha = (upperBound[i] + lowerBound[i]) / 2;
                for (ui j = laneId; j < total; j += warpSize){
                    ui neighborOffset = flowNetwork.neighborOffset2[fStart+j];
                    flowNetwork.capacity[neighborOffset] = alpha * k;


                }

        }

        }

    }
}
__global__ void getResult(deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp, deviceCliquesPointer finalCliqueData, ui *compCounter,  double * upperBound, double * lowerBound, double *densities, ui *checkResult, ui totalWarps, int totalComponents, ui k,ui t) {
  extern __shared__ char sharedMemory[];
  ui sizeOffset = 0;

  ui *size = (ui *)(sharedMemory + sizeOffset);
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int warpId = idx / warpSize;
  int laneId = idx % warpSize;

  for (ui i = warpId; i < totalComponents; i += totalWarps) {
    ui start = conComp.componentOffset[i];
    ui end = conComp.componentOffset[i + 1];
    ui total = end - start;
    ui cliqueStart = compCounter[i];
    ui cliqueEnd = compCounter[i + 1];
    ui totalCliques = cliqueEnd - cliqueStart;

    ui tFlow = total + totalCliques +2;

    //double bais = 1.0 / (tFlow * (tFlow - 1));


    if (checkResult[i]) {

        ui fStart = flowNetwork.offset[i];
        ui neighborOffset = flowNetwork.neighborOffset2[fStart+total+totalCliques];

        if(laneId==0){
            size[threadIdx.x/warpSize] =0;
        }
         __syncwarp();

        for(ui j = laneId; j < total; j+=warpSize){
            double residual = flowNetwork.capacity[neighborOffset+ j] - flowNetwork.flow[neighborOffset+j];
            if(residual>0){
                atomicAdd(&size[threadIdx.x/warpSize],1);
            }

        }
         __syncwarp();



        for(ui j = laneId; j < totalCliques; j +=warpSize){
            ui w =0;
            ui found = true;
            while(w<k){
                ui vertex = finalCliqueData.trie[w*t+j + cliqueStart ];
                double residual = flowNetwork.capacity[neighborOffset+ vertex] - flowNetwork.flow[neighborOffset+vertex];
                if(residual <= 0 ){
                    found = false;
                    break;
                }
            }

            if(found){
                atomicAdd(& densities[i],1);

            }

        }
         __syncwarp();
        if(laneId==0){
            //Density
            if(size[threadIdx.x/warpSize]==0){
              densities[i] = lowerBound[i];
            }else{
              densities[i]= densities[i]/size[threadIdx.x/warpSize];
            }
            //printf("warp %d size %d density %f \n",warpId,size[threadIdx.x/warpSize],densities[i]);



        }
         __syncwarp();


  }
}
}

__global__ void countCliques(deviceDAGpointer D, cliqueLevelDataPointer levelData,  ui *globalCounter, ui maxBitMask ,ui totalTasks, ui totalWarps){
    /* Find Total number of cliques in the graph by counting the valid neighbors.
       Each warp processes on partial clique. */

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(int i =warpId; i < totalTasks ; i+= totalWarps ){

        // candidate offset
        int start = levelData.offset[i];
        int totalCandidates = levelData.offset[i+1]- start;

        int count = 0;

        for(int j = laneId; j <totalCandidates; j+=warpSize){
            int degree = D.degree[levelData.candidates[start + j]];
            int numBitmasks = (degree + 31) / 32;
            for(int x =0; x <numBitmasks; x++){
                int neighBitMask = levelData.validNeighMask[(start+j)*maxBitMask + x];
                count += __popc(neighBitMask);

            }

        }

        __syncwarp();

         for (int offset = warpSize / 2; offset > 0; offset /= 2) {
                count += __shfl_down_sync(0xFFFFFFFF, count, offset);
        }

        if(laneId==0){
            atomicAdd(globalCounter,count);

        }



    }

}
