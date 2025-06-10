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
    /* Calculates the out degree of each vertex in DAG.
        One warp processes on vertex.*/
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    // Iter throught verticies. 
    for(ui i = warpId; i < n; i += totalWarps) {

        // Get the offset in orginal graph.
        ui start = G.offset[i];
        ui end = G.offset[i+1];

        // Total neighbors
        ui total = end - start;
        ui neigh;

        // stores total neighbors in DAG.
        int count = 0;

        // Each lanes iters through neighbors
        for(int j = laneId; j < total; j += warpSize) {
            neigh = G.neighbors[start + j];

            // If listing order is greater than vertex, edge is considered in DAG.
            if(listingOrder[i] < listingOrder[neigh]) {
                count++;
            }
        }

        // Reduce the count to get out degree in DAG.
        for (int offset = warpSize / 2; offset > 0; offset /= 2) {
            count += __shfl_down_sync(0xFFFFFFFF, count, offset);
        }

        // write in global memory
        if(laneId == 0) {
            D.degree[i] = count;
        }
    }
}

__global__ void generateNeighborDAG(deviceGraphPointers G, deviceDAGpointer D, ui *listingOrder, ui n, ui m, ui totalWarps) {

    /* Generate neighbor list of DAG based on the offset. 
       Each warp process on vertex at a time.
       Only threads (lanes) within a warp use atomicAdd to compete for writing to the same location.
       */

    // used to get the next write location in global memory.
    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui *counter = (ui *)(sharedMemory + sizeOffset);

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    // Iter through verticies
    for(ui i = warpId; i < n; i += totalWarps) {

        // Set counter to DAG offset of that vertex.
        if(laneId==0){
          counter[threadIdx.x / warpSize] = D.offset[i];
        }

        __syncwarp();

        // Offset in orginal graph
        ui start = G.offset[i];
        ui end = G.offset[i+1];

        // Total neghbors in orginal graph
        ui total = end - start;
        ui neigh;

        // Lanes iter through neighbors in orginal graph
        for(int j = laneId; j < total; j += warpSize) {
            neigh = G.neighbors[start + j];

            // If listing order of neighbor is greater than vertex, edge considered in DAG.
            if(listingOrder[i] < listingOrder[neigh]) {

                // Get next location and write the neighbor
                int loc = atomicAdd(&counter[threadIdx.x / warpSize], 1);
                D.neighbors[loc] = neigh;

            }
        }
      __syncwarp();
    }
}

__global__ void listIntialCliques(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui *label, ui k, ui n, ui m, ui psize, ui cpSize, ui maxBitMask, ui level, ui totalWarps) {
    /* 
        Generates initial partial cliques, their candidate extensions, and valid neighbor masks.
        The data is stored in virtual partitions, with each warp writing to a separate partition, to reduce memory contention.
        One warp processes one vertex and then moves to next.
    */
    
    // used to get next write location for a candidate in virtual partition.
    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui *counter = (ui *)(sharedMemory + sizeOffset);

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    // Virtual Partition offsets
    int cliquePartition = warpId * psize;
    int offsetPartition = warpId * (psize / (k-1) + 1);
    int candidatePartition = warpId * cpSize;
    int maskPartition = warpId * cpSize * maxBitMask;


    // Each Warp processes on vertex
    for(int i = warpId; i < n; i += totalWarps) {

        ui vertex = i;

        // DAG offset
        ui neighOffset = D.offset[vertex];


        if(laneId == 0) {
            counter[threadIdx.x / warpSize] = 0;
        }

        __syncwarp();


        int candidateOffset = candidatePartition + levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]];

        // Iter through neighbors of vertex and add to candidate in label = k.
        for(int j = laneId; j < D.degree[vertex]; j += warpSize) {
            ui neigh = D.neighbors[neighOffset + j];
            
            // Check and update the label. used to avoid duplication
            ui old_val = atomicCAS(&label[warpId*n + neigh], k, k-1);

            // if old label is k, i.e. neighbor yet not considered. write the neighbor in candidate array. 
            if(old_val==k) {

                ui loc = atomicAdd(&counter[threadIdx.x / warpSize], 1);
                levelData.candidatesPartition[candidateOffset + loc] = neigh;
            }
        }
        __syncwarp();

        // If atleast one valid candidate.
        if(laneId == 0 && counter[threadIdx.x / warpSize] > 0) {

            // Add vertex as a partial clique.
            levelData.partialCliquesPartition[cliquePartition + levelData.count[warpId + 1] * (k-1) + level] = vertex;
            
            // Increament number of partial cliques.
            levelData.count[warpId + 1] += 1;

            // Write candidate offset
            levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]] =
                levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1] - 1] + counter[threadIdx.x / warpSize];
        }
         __syncwarp();

        int start = candidateOffset;

        // Iter through each candidate, one lane processes on candidate.
        for(int j = laneId; j < counter[threadIdx.x / warpSize]; j += warpSize) {

            int candidate = levelData.candidatesPartition[start + j];
            int neighOffset = D.offset[candidate];
            int degree = D.degree[candidate];

            // Number of bitmask required to store the valid neighbor bitmask
            int numBitmasks = (degree + 31) / 32;

            // Generate bitmasks
            for (int bitmaskIndex = 0; bitmaskIndex < numBitmasks; bitmaskIndex++) {
                ui bitmask = 0; // Initialize bitmask to 0

                // Iterate over the current chunk of 32 neighbors
                int startNeighbor = bitmaskIndex * 32;
                int endNeighbor = min(startNeighbor + 32, degree);

                // Generate Bit Mask, if neighbor is valid.
                for (int x = startNeighbor; x < endNeighbor; x++) {
                    if (label[warpId*n + D.neighbors[neighOffset + x]] == k - 1) {
                        bitmask |= (1 << (x - startNeighbor)); // Set the bit for valid neighbors


                    }
                }
                // Write bitmask to memory
                levelData.validNeighMaskPartition[maskPartition + (levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]-1]+j) * maxBitMask + bitmaskIndex] = bitmask;
            }
        }

        __syncwarp();

        // update the label for next vertex. 
        for(int x = laneId; x<n;x+=32){
          label[warpId*n + x] = k;
        }

       __syncwarp();
    }

}

__global__ void flushParitions(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui pSize, ui cpSize, ui k, ui maxBitMask, ui level, ui totalWarps){
    /* Copy the partial cliques, candidates and valid neighbor bit masks from virtual partitions
         to contigeous arrays.
         One warp processes partial cliques of on virtual partition. */
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    // Virtual partition offsets.
    int cliquePartition = warpId * pSize;
    int offsetPartition = warpId * (pSize / (k-1) + 1);
    int candidatePartition = warpId * cpSize;
    int maskPartition = warpId * cpSize * maxBitMask;

    // Total number of partial cliques in a partition.
    int totalTasks = levelData.count[warpId+1] - levelData.count[warpId];

    // Iter through partial cliques in a virtual partition one at a time.
    for(int iter = 0; iter < totalTasks; iter++){

        // candidate offset on a partial clique.
        int start = candidatePartition + levelData.offsetPartition[offsetPartition + iter];
        int end = candidatePartition + levelData.offsetPartition[offsetPartition + iter+ 1];

        // Total number of candidates.
        int total = end-start;

        // write offset of candidates in the contigeous array.
        int writeOffset = levelData.temp[warpId] + levelData.offsetPartition[offsetPartition + iter];

        // Each lane processes one candidate parallely
        for(int i = laneId; i < total; i+=warpSize){
           ui candidate = levelData.candidatesPartition[start + i];

           // write in contigeous memory.
           levelData.candidates[writeOffset+ i] = candidate;

           // Get total mask required to store the valid bitmasks.
            int totalMasks = (D.degree[candidate]+31)/32;

            // Write the valid bit masks in contigeous memory.
            for(int j =0; j < totalMasks; j++){
                levelData.validNeighMask[(writeOffset+i)*maxBitMask + j ] =
                levelData.validNeighMaskPartition[maskPartition + (levelData.offsetPartition[offsetPartition + iter] + i)*maxBitMask + j];
            }

        }

        // first level lanes, write the partial clique in contigeous memory. 
        if(laneId< level+1 ){
            levelData.partialCliques[levelData.count[warpId]*(k-1)+ iter*(k-1) + laneId] = levelData.partialCliquesPartition[cliquePartition + iter * (k-1) + laneId];
        }

        __syncwarp();

        // update the candidate offset in contigeous array.
        if(laneId==0){
            levelData.offset[levelData.count[warpId] + iter + 1] = levelData.temp[warpId]+levelData.offsetPartition[offsetPartition + iter+ 1];
        }
        __syncwarp();


    }

}

__global__ void listMidCliques(deviceDAGpointer D, cliqueLevelDataPointer levelData, ui *label, ui k,ui iterK, ui n, ui m,ui pSize, ui cpSize, ui maxBitMask,ui totalTasks, ui level, ui totalWarps){
    /* This kernel extends partial cliques by adding new vertices and updates the candidate set along with their valid neighbors.
       Each warp processes the one partial cliques at a time.
        */
    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui *counter = (ui * )(sharedMemory + sizeOffset);

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    // Virtual partition offsets
    int cliquePartition  = warpId*pSize;
    int offsetPartition = warpId*(pSize/(k-1)+1);
    int candidatePartition = warpId*cpSize;
    int maskPartition = warpId*cpSize*maxBitMask;

    // Iter through partial cliques.
    for(int i =warpId; i < totalTasks ; i+= totalWarps ){

        // Candidate offset
        int start = levelData.offset[i];

        int totalCandidates = levelData.offset[i+1]- start;

        // Process one canddate at a time
        for(int iter = 0; iter <totalCandidates; iter ++){

            int candidate = levelData.candidates[start + iter];

            if(laneId==0){
                counter[threadIdx.x/warpSize] = 0;
            }
            __syncwarp();

            // Candidate degree in DAG
            int degree = D.degree[candidate];
            int neighOffset = D.offset[candidate];

            // write offset in virtual partition. 
            int writeOffset = candidatePartition + levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]];
            
            // Each lane processes neighbors.
            for(int j = laneId; j< degree; j+= warpSize ){

                int iterBitMask = j/warpSize;
                int bitPos = j%32;
                int neighBitMask = levelData.validNeighMask[(start+iter)*maxBitMask + iterBitMask];
                ui neigh = D.neighbors[neighOffset + j];

                // if neighbor is valid.
                if(neighBitMask & (1 << bitPos )){

                     // compare and update the label.
                     ui old_val = atomicCAS(&label[warpId*n + neigh], iterK, iterK-1);

                     // if visited first time.
                     if(old_val== iterK){

                        // write vertex as a new candidate
                        ui loc = atomicAdd( &counter[threadIdx.x/warpSize], 1);
                        levelData.candidatesPartition[writeOffset + loc] = neigh;

                     }
                }

            }
            __syncwarp();
            if(laneId == 0 && counter[threadIdx.x/warpSize] > 0){

                // current candidate is added to the partial clique.
                levelData.partialCliquesPartition[cliquePartition + levelData.count[warpId+1] * (k-1) + level ] = candidate;

                // copy the rest of the partial clique.
                for(int l =0; l<level; l++){
                  levelData.partialCliquesPartition[cliquePartition + levelData.count[warpId+1] * (k-1) + l ] = levelData.partialCliques[i*(k-1)+l];
                }

                // Increament pc count.
                levelData.count[warpId+1] +=1;

                // Update pc offset.
                levelData.offsetPartition[offsetPartition + levelData.count[warpId+1]] =
                    levelData.offsetPartition[offsetPartition + levelData.count[warpId+1] - 1] +counter[threadIdx.x/warpSize];
            }

            __syncwarp();

            int start = writeOffset;

            // iter through candidates.
            for(int j = laneId; j < counter[threadIdx.x / warpSize]; j += warpSize) {

                int cand = levelData.candidatesPartition[start + j];
                int neighOffset = D.offset[cand];
                int degree = D.degree[cand];

                // total bit masks required
                int numBitmasks = (degree + 31) / 32;

                // 
                for (int bitmaskIndex = 0; bitmaskIndex < numBitmasks; bitmaskIndex++) {
                    ui bitmask = 0; // Initialize bitmask to 0

                    // Iterate over the current chunk of 32 neighbors
                    int startNeighbor = bitmaskIndex * 32;
                    int endNeighbor = min(startNeighbor + 32, degree);
                    for (int x = startNeighbor; x < endNeighbor; x++) {


                        if (label[warpId*n + D.neighbors[neighOffset + x]] == (iterK - 1)) {
                            bitmask |= (1 << (x - startNeighbor)); // Set the bit for valid neighbors


                        }
                    }

                    // write bit mask in global memory
                    levelData.validNeighMaskPartition[maskPartition + (levelData.offsetPartition[offsetPartition + levelData.count[warpId + 1]-1]+j) * maxBitMask + bitmaskIndex] = bitmask;
                }
            }

            __syncwarp();

            // update labels for next pc.
            for(int x = laneId; x<n;x+=32){
              label[warpId*n + x] = iterK;
            }


           __syncwarp();
        }

    }

}

__global__ void writeFinalCliques(deviceGraphPointers G, deviceDAGpointer D, cliqueLevelDataPointer levelData, deviceCliquesPointer cliqueData, ui *globalCounter,ui k,ui iterK, ui n, ui m,ui pSize, ui cpSize, ui maxBitMask,ui trieSize,ui totalTasks, ui level, ui totalWarps){
    /* Write final k-cliques from k-2 partial cliques.
       Each warp processes on partial clique. */
      
    
    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;
    ui *counter = (ui * )(sharedMemory + sizeOffset);

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    for(int i =warpId; i < totalTasks ; i+= totalWarps ){

        // candidate offset
        int start = levelData.offset[i];
        int totalCandidates = levelData.offset[i+1]- start;

        // iter through candidates sequentially. 
        for(int iter = 0; iter <totalCandidates; iter ++){

            int candidate = levelData.candidates[start + iter];
            if(laneId==0){
                counter[threadIdx.x/warpSize]=0;
            }
            __syncwarp();
            int degree = D.degree[candidate];
            int neighOffset = D.offset[candidate];

            // Iter through neighbors of candidate
            for(int j = laneId; j< degree; j+= warpSize ){

                int iterBitMask = j/warpSize;
                int bitPos = j%32;
                int neighBitMask = levelData.validNeighMask[(start+iter)*maxBitMask + iterBitMask];

                // if neighbor is valid.
                if(neighBitMask & (1 << bitPos )){

                    ui neigh = D.neighbors[neighOffset + j];

                    // copy the partial clique to final clique data.
                    ui loc = atomicAdd(globalCounter,1);
                    for(int ind =0; ind < k-2; ind++){
                        cliqueData.trie[trieSize * ind + loc] = levelData.partialCliques[(i)*(k-1) + ind];

                    }
                    
                    //increament Clique degree of the all verticies in partial clique.
                    atomicAdd(&counter[threadIdx.x/warpSize],1);

                    // add candidate to clique
                    cliqueData.trie[trieSize * (k-2) + loc]  = candidate;
                    // add neighbor to clique
                    cliqueData.trie[trieSize * (k-1) + loc] = neigh;
                    // set status of clique to -1 i.e. valid clique
                    cliqueData.status[loc]= -1;
                    // increase clique degree of neigh and candidate
                    atomicAdd(&G.cliqueDegree[neigh],1);
                    atomicAdd(&G.cliqueDegree[candidate],1);

                }


            }
            __syncwarp();

            // update clique degree of all verticies in partial clique.
            for(int j = laneId; j< k-2 ; j+= warpSize ){
                int pClique = levelData.partialCliques[i*(k-1) + j];
                atomicAdd(&G.cliqueDegree[pClique],counter[threadIdx.x/warpSize]);
            }
             __syncwarp();

        }
    }

}


/*__global__ void sortTrieData(deviceGraphPointers G, deviceCliquesPointer cliqueData, ui totalCliques, ui t, ui k, ui totalThreads){
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

}*/


__global__ void selectNodes(deviceGraphPointers G, ui *bufTails,ui *glBuffers, ui glBufferSize, ui n, ui level){
    __shared__ ui *glBuffer;
    __shared__ ui bufTail;

    // get the glbuffer for this warp.
    if(threadIdx.x == 0){
        bufTail = 0;
        glBuffer = glBuffers + blockIdx.x*glBufferSize;
    }
    __syncthreads();

    ui idx = blockIdx.x * blockDim.x + threadIdx.x;

    // iter through the verticies
    for(ui i = idx ;i<n; i+=BLK_DIM){
      ui v = i;
      // if core value is equal to level.
      if(G.cliqueCore[v] == level){

        // add to its coresponding glbuffer and increament bufftail.
        ui loc = atomicAdd(&bufTail, 1);
        glBuffer[loc] = v;


        }

    }


    __syncthreads();

    if(threadIdx.x == 0)
    {
        //printf("idx %d block %d bufTail %d \n", idx,blockIdx.x,bufTail);
        bufTails [blockIdx.x] = bufTail;
    }

}

__global__ void processNodesByWarp(deviceGraphPointers G,deviceCliquesPointer cliqueData, ui *bufTails,ui *glBuffers, ui *globalCount, ui glBufferSize, ui n, ui level, ui k, ui t, ui tt){
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
            if(regTail < base )
            base = regTail;
        }

        //vertex to be removed
        ui v = glBuffer[i];


        __syncwarp();

        // warp iters through the clique data.
        for(ui j =laneId; j<tt; j+=warpSize){

            // if valid clique and not removed yet
            if(cliqueData.status[j] == -1){

                // flag to check if vertex found in the clique
                bool found = false;
                // stores the index at which the vertex was found in clique (0,k-1)
                ui w =0;

                // iter through verticies of clique sequentially. 
                while(w<k){

                    if(cliqueData.trie[w*t+j] == v){

                    found = true;
                    break;
                    }
                    w++;
                }

                if(found){
                    // iter throught the clique verticies
                    for(ui x =0 ;x<k;x++){

                        // continue it clique vertex is same as vertex
                        if(x==w) continue;
                        // clique vertex
                        ui u = cliqueData.trie[x*t+j];
                        // decreament its core value by 1
                        int a = atomicSub(&G.cliqueCore[u],1);
                        
                        // if core value is less than level, update to level.
                        if (a <= level) {
                            atomicMax(&G.cliqueCore[u], level);
                        }
                        // if core value is level, add to glbuffer so can be removed in this level.
                        if(a == level+1){
                            ui loc = atomicAdd(&bufTail, 1);
                            glBuffer[loc] = u;

                        }
                    }

                    // set status of the clique to current core level.
                    cliqueData.status[j] = level;


                }
            }
        }

        __syncwarp();
        if(laneId == 0 && bufTail>0){
            // increament global count as one vertex was removed
            atomicAdd(globalCount, 1);
        }
    }
}


// TODO: need update
__global__ void processNodesByBlock(deviceGraphPointers G,deviceCliquesPointer cliqueData, ui *bufTails,ui *glBuffers, ui *globalCount, ui glBufferSize, ui n, ui level, ui k, ui t, ui tt){
    /*removes the verties to get their core value.
      block processes the verticies in its virtual partition */
    
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

            if( (v == cliqueData.trie[j]) && (cliqueData.status[j] == -1 )){
                for(ui x =1;x<k;x++){
                    ui u = cliqueData.trie[x*t+j];
                    int a = atomicSub(&G.cliqueCore[u], 1);
                    if(a == level+1){
                        ui loc = atomicAdd(&bufTail, 1);
                        glBuffer[loc] = u;

                    }
                    if(a <= level){
                        atomicMax(&G.cliqueCore[u], level);
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

__global__ void generateDensestCore(deviceGraphPointers G, densestCorePointer densestCore,ui *globalCount, ui n, ui density, ui totalWarps){
    /*Generates densest core, remaps it and calculates the new offset.
      Each warp processes on vertex at a time.*/
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    // iter through verticies
    for(ui i = warpId; i < n; i += totalWarps){
        // if core value is greater or equal to core value of densest core.
        if(G.cliqueCore[i]>= density){
            ui loc;
            // add vertex to densest core. 
            if(laneId==0){
                loc = atomicAdd(globalCount,1);
                densestCore.mapping[loc] = i;
            }
            // so loc can be used by other lanes
            loc = __shfl_sync(0xFFFFFFFF, loc, 0, 32);

            // Neighbor details
            ui start = G.offset[i];
            ui end = G.offset[i+1];
            ui total = end - start;
            ui neigh;
            int count = 0;
            // count neighbors who density is also greater
            for(int j = laneId; j < total; j += warpSize) {
                neigh = G.neighbors[start + j];
                if(G.cliqueCore[neigh] >= density) {
                    count++;
                }
            }
            // reduce to get total
            for (int offset = warpSize / 2; offset > 0; offset /= 2) {
                count += __shfl_down_sync(0xFFFFFFFF, count, offset);
            }
            __syncwarp();

            if(laneId == 0) {
                // update offset
                densestCore.offset[loc+1] = count;
            }
            __syncwarp();

        }
    }
}

__global__ void generateNeighborDensestCore(deviceGraphPointers G, densestCorePointer densestCore, ui density, ui totalWarps) {
    /* Generates neighbors of each vertex in densest core.
      each warp processes a vertex in densest core*/
    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui *counter = (ui *)(sharedMemory + sizeOffset);

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    // iter through verticies of densest core
    for(ui i = warpId; i < (*densestCore.n); i += totalWarps) {
        if(laneId==0){
          // set counter to new offset.
          counter[threadIdx.x / warpSize] = densestCore.offset[i];
        }
        __syncwarp();
        ui vertex = densestCore.mapping[i];

        // Neighbors in orginal graph
        ui start = G.offset[vertex];
        ui end = G.offset[vertex+1];
        ui total = end - start;
        ui neigh;
        
        for(int j = laneId; j < total; j += warpSize) {
            neigh = G.neighbors[start + j];
            // if clique core is greater
            if(G.cliqueCore[neigh] >= density) {
                int loc = atomicAdd(&counter[threadIdx.x / warpSize], 1);
                // write remapped neighbor in new neighbor list.
                densestCore.neighbors[loc] = densestCore.reverseMap[neigh];

            }
        }
      __syncwarp();
    }
}

__global__ void pruneEdges(densestCorePointer densestCore, deviceCliquesPointer cliqueData, ui *pruneStatus,ui t, ui tt, ui k, ui level ){
    /* Removes edges from the densest core that are not part of any clique.
     Each warp processes on clique. */
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    // iter through all the cliques.
    for(ui i = warpId; i<tt; i+=TOTAL_WARPS){

        // if clique is part of the densest core.
        if(cliqueData.status[i] >= level){

            // iter through each edge sequentially 
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

                        // find the edge is the neighbor list and update its status
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
    /*calculates new degree of the densest core after the edge pruning
      each warp process one vertex.*/
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    // Iter through each vertex
    for(ui i = warpId; i < n; i += totalWarps) {

        // Neighbors in orginal densest core.
        ui start = densestCore.offset[i];
        ui end = densestCore.offset[i+1];
        ui total = end - start;
        int count = 0;

        // increament count if not prunned
        for(int j = laneId; j < total; j += warpSize) {
            if(!pruneStatus[start + j]) {
                count++;
            }
        }

        // Reduce to get total new neighbors
        for (int offset = warpSize / 2; offset > 0; offset /= 2) {
            count += __shfl_down_sync(0xFFFFFFFF, count, offset);
        }

        if(laneId == 0) {
            // write new offset
            newOffset[i+1] = count;
        }
    }
}

__global__ void generateNeighborAfterPrune(densestCorePointer densestCore ,ui *pruneStatus, ui *newOffset, ui *newNeighbors,ui n, ui m, ui totalWarps) {
    /* Generate new neighbor list after pruning. 
       Each Warp processes one vertex*/
    extern __shared__ char sharedMemory[];
    ui sizeOffset = 0;

    ui *counter = (ui *)(sharedMemory + sizeOffset);

    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;

    // iter through densest core 
    for(ui i = warpId; i < n; i += totalWarps) {
        if(laneId==0){
          counter[threadIdx.x / warpSize] = newOffset[i];
        }
        __syncwarp();
        ui start = densestCore.offset[i];
        ui end = densestCore.offset[i+1];
        ui total = end - start;
        ui neigh;

        // write neighbors is not pruned
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
    /* Uses convergence to decompose densest core into its connected components.
      one warp process one vertex at a time.*/
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int warpId = idx / warpSize;
    int laneId = idx % warpSize;
    bool threadChanged = false;

    // iter through all verticies of the densest core.
    for(ui i = warpId; i < n; i += totalWarps) {

        // current component id.
        ui currentComp = conComp.components[i];

        // Neighbors 
        ui start = prunedNeighbors.newOffset[i];
        ui end = prunedNeighbors.newOffset[i+1];
        ui total = end - start;

        ui minNeighComp = currentComp;

        // iter through neighbors
        for (ui j = laneId; j < total; j += warpSize) {
            // Find minimum of component id amoung neighbor and vertex
            ui neighComp = conComp.components[prunedNeighbors.newNeighbors[start+j]];
            minNeighComp = min(minNeighComp, neighComp);
        }

        // Reduce to find minimum among all neighbors
        for (int offset = warpSize/2; offset > 0; offset /= 2) {
            ui temp = __shfl_down_sync(0xFFFFFFFF, minNeighComp, offset);
            minNeighComp = min(minNeighComp, temp);
        }

        if (laneId == 0) {
            // update the component id of vertex to mimimum component id.
            if ( minNeighComp < currentComp) {
                conComp.components[i] = minNeighComp;
                // set flag to true
                threadChanged = true;

            }
        }

        __syncwarp();
    }

    // if any lane in warp changed the component id, set flag to true.
    bool warpChanged = __any_sync(0xFFFFFFFF, threadChanged);

    // increament changed to indicate that component id was changed
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
            comp = min(comp,conComp.components[conComp.reverseMapping[vertex]]);

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


__global__ void createFlowNetworkOffset(deviceGraphPointers G, deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp, densestCorePointer densestCore, deviceCliquesPointer finalCliqueData, ui *compCounter,double *upperBound , ui totalWarps, ui totalComponents, ui k, double lb, ui t){



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

__global__ void createFlowNetwork(deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp, densestCorePointer densestCore, deviceCliquesPointer finalCliqueData, ui *compCounter,double *upperBound , ui totalWarps, ui totalComponents, ui k, double lb, ui t){



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

            double alpha = upperBound[i]/lb;

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

__global__ void pushRelabel(deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp, densestCorePointer densestCore, deviceCliquesPointer finalCliqueData, ui * compCounter, double * upperBound, double * lowerBound, ui * activeNodes, ui * componenetsLeft, ui totalWarps, int totalComponents, ui k,ui t, ui partitionSize) {
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

      int maxIterations = 4;

      //printf("warpid %d laneid %d\n ", warpId,laneId);

      //Push or Relabel until converge
      for (int iter = 0; iter < maxIterations; iter++) {
        if (laneId == 0) {
          counter[threadIdx.x / warpSize] = 0;
        }
        __syncwarp();

        //check for nodes with excess
        for (ui j = laneId; j < tFlow; j += warpSize) {
          if (j != s && j != t && flowNetwork.excess[fStart + j] > 0) {
            int pos = atomicAdd( & counter[threadIdx.x / warpSize], 1);
            if (pos < partitionSize) {
              activeNodes[i * partitionSize + pos] = j;
            }
           // printf("warp id %d lane id %d pos %d u %d  iter %d \n", i, j, pos, j, iter);

          }
        }
        __syncwarp();
        bool should_break = (counter[threadIdx.x / warpSize] == 0);
        __syncwarp();  // Sync before any thread exits

        if (should_break) break;

        bool pushed = false;

        for (ui j = laneId; j < counter[threadIdx.x / warpSize]; j += warpSize) {
          ui vertex = activeNodes[i * partitionSize + j];

          ui nStart = flowNetwork.neighborOffset2[fStart + vertex];
          ui nEnd = flowNetwork.neighborOffset2[fStart + vertex + 1];
          //printf("warpid %d lane id %d read ver %d excess %f height %d offset %d nstart %d end %d \n", i, j, vertex, flowNetwork.excess[fStart + vertex], flowNetwork.height[fStart + vertex], fStart, nStart, nEnd);

          //Check neighbors to send acess to.
          for (ui x = nStart; x < nEnd; x++) {
            //
            /*if (flowNetwork.excess[fStart + vertex] == 0) {
              break;
            }*/
            ui neigh = flowNetwork.Edges[x];
            double residual = flowNetwork.capacity[x] - flowNetwork.flow[x];


            __syncwarp();

            // If neighbor has capacity
            if ((flowNetwork.height[fStart + vertex] == flowNetwork.height[fStart + neigh] + 1) && residual > 0) {

              //printf("warp Id %d laneid %d vertex %d nstart %d end %d neigh %d resi %f cap %f flow %f loc %d hv %d hu %d ev %f en %f \n", i, j, vertex, nStart, nEnd, neigh, residual, flowNetwork.capacity[x], flowNetwork.flow[x], x, flowNetwork.height[fStart + vertex], flowNetwork.height[fStart + neigh], flowNetwork.excess[fStart + vertex], flowNetwork.excess[fStart + neigh]);

              double delta = min(flowNetwork.excess[fStart + vertex], residual);
              if (delta > 0) {

                //forward flow vertex to neigh
                atomicAdd( & flowNetwork.flow[x], delta);

                //Backward Flow neigh to vertex
                ui stemp = flowNetwork.neighborOffset2[fStart + neigh];
                ui etemp = flowNetwork.neighborOffset2[fStart + neigh + 1];
                for (ui ind = stemp; ind < etemp; ind++) {
                  if (flowNetwork.Edges[ind] == vertex) {
                    atomicAdd( & flowNetwork.flow[ind], -delta);
                    break;
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

          __syncwarp();

          // Relabel
          if (!pushed && flowNetwork.excess[fStart + vertex] > 0) {
            ui minHeight = UINT_MAX;
            //printf("Relabel  WarpId %d laneid %d vertex %d excess %f height %d nstart %d end %d \n", i, laneId, vertex, flowNetwork.excess[fStart + vertex], flowNetwork.height[fStart + vertex], nStart, nEnd);
            for (ui x = nStart; x < nEnd; x++) {
              ui neigh = flowNetwork.Edges[x];
              double residual = flowNetwork.capacity[x] - flowNetwork.flow[x];
              //printf("warp %d laneid %d j %d v %d res %f negh %d height %d nstart %d \n", i, laneId, j, vertex, residual, neigh, flowNetwork.height[fStart + neigh], nStart);

              if (residual > 0.0) {
                minHeight = (flowNetwork.height[fStart + neigh] < minHeight) ? flowNetwork.height[fStart + neigh] : minHeight;
              }

            }
            if (minHeight != INF) {
              flowNetwork.height[fStart + vertex] = minHeight + 1;

            }

            //printf("after Relabel  WarpId %d laneid %d vertex %d excess %f height %d minHeight %d \n", i, laneId, vertex, flowNetwork.excess[fStart + vertex], flowNetwork.height[fStart + vertex], minHeight);

          }
          __syncwarp();
        }

      }

      __syncwarp();


      //printf("second warpid %d laneid %d\n ", warpId,laneId);

      if(laneId==0){
        //getRes[threadIdx.x / warpSize] = 0;
        double alpha = (upperBound[i] + lowerBound[i]) / 2;
        if (flowNetwork.excess[fStart + t] == (double) totalCliques * k) {
            upperBound[i] = alpha;
        } else {
            lowerBound[i] = alpha;

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
__global__ void getResult(deviceFlowNetworkPointers flowNetwork, deviceComponentPointers conComp, deviceCliquesPointer finalCliqueData, ui *compCounter,  double * upperBound, double * lowerBound, double *densities,ui totalWarps, int totalComponents, ui k,ui t) {
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

    double bais = 1.0 / (tFlow * (tFlow - 1));


    if (((upperBound[i] - lowerBound[i]) < bais)&&(upperBound[i]!=0)&&(lowerBound[i]!=0)) {

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
        

        }
         __syncwarp();


  }
}
}


