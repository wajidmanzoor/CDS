#include "../inc/helpers.cuh"
#include "../utils/cuda_utils.cuh"

#define ASSERT_DEBUG(cond, fmt, ...)                                           \
  if (!(cond)) {                                                               \
    printf("ASSERT FAILED at %s:%d: " fmt, __FILE__, __LINE__, ##__VA_ARGS__); \
    assert(false);                                                             \
  }

__global__ void generateDegreeDAG(deviceGraphPointers G, deviceDAGpointer D,
                                  ui *listingOrder, ui n, ui m, ui totalWarps) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int warpId = idx / warpSize;
  int laneId = idx % warpSize;

  for (ui i = warpId; i < n; i += totalWarps) {
    ui start = G.offset[i];
    ui end = G.offset[i + 1];
    ui total = end - start;
    ui neigh;
    int count = 0;
    for (int j = laneId; j < total; j += warpSize) {
      neigh = G.neighbors[start + j];
      if (listingOrder[i] < listingOrder[neigh]) {
        count++;
      }
    }

    for (int offset = warpSize / 2; offset > 0; offset /= 2) {
      count += __shfl_down_sync(0xFFFFFFFF, count, offset);
    }

    if (laneId == 0) {
      D.degree[i] = count;
    }
  }
}

__global__ void generateNeighborDAG(deviceGraphPointers G, deviceDAGpointer D,
                                    ui *listingOrder, ui n, ui m,
                                    ui totalWarps) {

  extern __shared__ char sharedMemory[];
  ui sizeOffset = 0;

  ui *counter = (ui *)(sharedMemory + sizeOffset);

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int warpId = idx / warpSize;
  int laneId = idx % warpSize;

  for (ui i = warpId; i < n; i += totalWarps) {
    if (laneId == 0) {
      counter[threadIdx.x / warpSize] = D.offset[i];
    }
    __syncwarp();
    ui start = G.offset[i];
    ui end = G.offset[i + 1];
    ui total = end - start;
    ui neigh;
    for (int j = laneId; j < total; j += warpSize) {
      neigh = G.neighbors[start + j];

      if (listingOrder[i] < listingOrder[neigh]) {
        int loc = atomicAdd(&counter[threadIdx.x / warpSize], 1);
        D.neighbors[loc] = neigh;
      }
    }
    __syncwarp();
  }
}

__global__ void listInitialCliquesBaseline(deviceDAGpointer D,
                                           cliqueLevelDataBaseline levelData,
                                           ui *label, ui k, ui n, ui maxBitMask,
                                           ui totalWarps) {
  extern __shared__ char sharedMemory[];
  ui *counter = (ui *)sharedMemory;

  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  int warpId = idx / warpSize;
  int laneId = idx % warpSize;
  ui oneLabelSize = (n + 31) / 32;
  size_t labelWordOffset = (size_t)warpId * (size_t)oneLabelSize;
  ui *warpLabel = (ui *)label + labelWordOffset;

  for (ui i = warpId; i < n; i += totalWarps) {

    ui vertex = i, vOff = D.offset[vertex];
    ui localTaskCount = -1;
    if (laneId == 0) {
      counter[threadIdx.x / warpSize] = 0;

      // Spin until we acquire the lock
      while (atomicCAS(&lock, 0, 1) != 0) {
        // busy wait
      }

      // Critical section
      localTaskCount = atomicAdd(&levelData.taskCount, 1);

      levelData.partialCliques[localTaskCount * (k - 1) + 0] = vertex;
      levelData.offset[localTaskCount + 1] =
          (localTaskCount == 0)
              ? D.degree[vertex]
              : levelData.offset[localTaskCount] + D.degree[vertex];

      // Release lock
      __threadfence(); // ensure writes are visible
      atomicExch(&lock, 0);
    }

    // Broadcast taskCount to all threads in the warp
    taskCount = __shfl_sync(0xFFFFFFFF, localTaskCount, 0);
    ui cOff = levelData.offset[taskCount];
    for (ui j = laneId; j < D.degree[vertex]; j += warpSize) {
      ui neigh = D.neighbors[vOff + j];
      ui wordIdx = neigh / 32; // Which ui word (0, 1, 2, ...)
      ui bitPos = neigh % 32;  // Which bit in that word (0-31)
      ui mask = 1U << bitPos;

      ui old = atomicOr(&warpLabel[wordIdx], mask);
      /*printf(" warp id %d i %d v %d j %d neigh %d deg %d vertex offset %d "
             "wordIdx %d bitPos %d old %d  \n",
             warpId, i, vertex, j, neigh, D.degree[vertex], vOff, wordIdx,
             bitPos, old);*/
      if (!(old & mask)) {
        ui loc = atomicAdd(&counter[threadIdx.x / warpSize], 1);
        levelData.candidates[cOff + loc] = neigh;
      }
    }
    __syncwarp();
    ui counterValue = counter[threadIdx.x / warpSize];
    for (ui j = laneId; j < counterValue; j += warpSize) {
      ui cand = levelData.candidates[cOff + j];
      ui dOff = D.offset[cand], deg = D.degree[cand];
      ui chunks = (deg + 31) / 32;

      for (ui m = 0; m < chunks; m++) {
        ui bitmask = 0, from = m * 32;
        ui to = min(from + 32, deg);
        for (ui x = from; x < to; x++) {
          ui nb = D.neighbors[dOff + x];
          // 32-bit atomic bit check
          ui wordIdx = nb / 32;
          ui bitPos = nb % 32;
          ui mask = 1U << bitPos;

          if ((warpLabel[wordIdx] & mask) != 0)
            bitmask |= 1 << (x - from);
        }

        // Use size_t for large index calculations
        size_t maskIdx = (cOff + j) * (size_t)maxBitMask + m;
        levelData.validNeighMask[maskIdx] = bitmask;
        /*printf(" warp id %d i %d v %d valid location %zu mask %d cand %d \n",
               warpId, i, vertex, maskIdx, bitmask, cand);*/
      }
    }
    __syncwarp();

    for (ui x = laneId; x < n; x += warpSize) {
      ui wordIdx = x / 32;
      ui bitPos = x % 32;
      ui mask = ~(1U << bitPos); // Inverted mask to clear the bit
      atomicAnd(&warpLabel[wordIdx], mask);
    }

    __syncwarp();
  }
}

__global__ void listMidCliquesBaseline(deviceDAGpointer D,
                                       cliqueLevelDataBaseline levelDataRead,
                                       cliqueLevelDataBaseline levelDataWrite,
                                       ui *label, ui k, ui n, ui maxBitMask,
                                       ui level, ui totalWarps) {
  extern __shared__ char sharedMemory[];
  ui *counter = (ui *)sharedMemory;

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int warpId = idx / warpSize;
  int laneId = idx % warpSize;

  ui oneLabelSize = (n + 31) / 32;
  ui *warpLabel = label + (size_t)warpId * oneLabelSize;
  ui totalTasks = levelDataRead.taskCount;

  // Iterate over tasks from previous level
  for (ui t = warpId; t < totalTasks; t += totalWarps) {

    ui start = levelDataRead.offset[t];
    ui end = levelDataRead.offset[t + 1];
    ui numCandidates = end - start;

    // Expand each candidate of this partial clique
    for (ui it = 0; it < numCandidates; it++) {

      ui pivot = levelDataRead.candidates[start + it];
      ui deg = D.degree[pivot];
      ui pOff = D.offset[pivot];
      ui localCount = 0;

      for (ui j = laneId; j < deg; j += warpSize) {
        ui maskIdx = j / 32;
        ui bitPos = j % 32;

        ui validMask =
            levelDataRead.validNeighMask[(start + it) * maxBitMask + maskIdx];

        if (validMask & (1U << bitPos)) {
          ui nb = D.neighbors[pOff + j];

          ui wordIdx = nb / 32;
          ui bit = nb % 32;
          ui mask = 1U << bit;

          if (!(warpLabel[wordIdx] & mask)) {
            atomicAdd(&localCount, 1);
          }
        }
      }

      // Reduce count across warp
      for (int offset = 16; offset > 0; offset >>= 1)
        localCount += __shfl_down_sync(0xffffffff, localCount, offset);

      ui writeTask = 0;
      ui writeOffset = 0;

      // --------------------------------
      // Phase 2: reserve output space
      // --------------------------------
      if (laneId == 0 && localCount > 0) {

        while (atomicCAS(&lock, 0, 1) != 0) {
        }

        writeTask = atomicAdd(&levelDataWrite.taskCount, 1);

        levelDataWrite.offset[writeTask + 1] =
            levelDataWrite.offset[writeTask] + localCount;

        atomicExch(&lock, 0);
      }

      writeTask = __shfl_sync(0xffffffff, writeTask, 0);
      writeOffset = levelDataWrite.offset[writeTask];

      if (localCount == 0)
        continue;

      // --------------------------------
      // Phase 3: write candidates
      // --------------------------------

      ui counter = 0;

      for (ui j = laneId; j < deg; j += warpSize) {
        ui maskIdx = j / 32;
        ui bitPos = j % 32;

        ui validMask =
            levelDataRead.validNeighMask[(start + it) * maxBitMask + maskIdx];

        if (validMask & (1U << bitPos)) {
          ui nb = D.neighbors[pOff + j];

          ui wordIdx = nb / 32;
          ui bit = nb % 32;
          ui mask = 1U << bit;

          ui old = atomicOr(&warpLabel[wordIdx], mask);
          if (!(old & mask)) {
            ui pos = atomicAdd(&counter, 1);
            levelDataWrite.candidates[writeOffset + pos] = nb;
          }
        }
      }

      __syncwarp();

      // --------------------------------
      // Write partial clique
      // --------------------------------
      if (laneId == 0) {
        size_t base = (size_t)writeTask * (k - 1);

        for (ui l = 0; l < level; l++) {
          levelDataWrite.partialCliques[base + l] =
              levelDataRead.partialCliques[(size_t)t * (k - 1) + l];
        }
        levelDataWrite.partialCliques[base + level] = pivot;
      }

      __syncwarp();

      // --------------------------------
      // Build valid neighbor mask
      // --------------------------------
      for (ui j = laneId; j < localCount; j += warpSize) {
        ui cand = levelDataWrite.candidates[writeOffset + j];
        ui dOff = D.offset[cand];
        ui dDeg = D.degree[cand];

        ui chunks = (dDeg + 31) / 32;

        for (ui m = 0; m < chunks; m++) {
          ui bitmask = 0;
          ui from = m * 32;
          ui to = min(from + 32, dDeg);

          for (ui x = from; x < to; x++) {
            ui nb = D.neighbors[dOff + x];
            ui w = nb / 32;
            ui b = nb % 32;
            if (warpLabel[w] & (1U << b))
              bitmask |= (1U << (x - from));
          }

          levelDataWrite.validNeighMask[(writeOffset + j) * maxBitMask + m] =
              bitmask;
        }
      }

      __syncwarp();

      // --------------------------------
      // Clear label
      // --------------------------------
      for (ui x = laneId; x < n; x += warpSize) {
        ui wordIdx = x / 32;
        ui bitPos = x % 32;
        ui mask = ~(1U << bitPos); // Inverted mask to clear the bit
        atomicAnd(&warpLabel[wordIdx], mask);
      }

      __syncwarp();
    }
  }
}

__global__ void writeFinalCliquesBaseline(cliqueLevelDataBaseline levelData,
                                          deviceDAGpointer D,
                                          deviceCliquesPointer C,
                                          ui *globalCounter, ui k) {

  int warpId = (blockIdx.x * blockDim.x + threadIdx.x) / warpSize;
  int laneId = threadIdx.x % warpSize;

  extern __shared__ char sharedMemory[];
  ui sizeOffset = 0;
  ui *counter = (ui *)(sharedMemory + sizeOffset);

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int warpId = idx / warpSize;
  int laneId = idx % warpSize;
  ui totalTasks = levelData.taskCount;
  for (int i = warpId; i < totalTasks; i += totalWarps) {

    int start = levelData.offset[i];
    int totalCandidates = levelData.offset[i + 1] - start;
    for (int iter = 0; iter < totalCandidates; iter++) {
      int candidate = levelData.candidates[start + iter];
      if (laneId == 0) {
        counter[threadIdx.x / warpSize] = 0;
      }
      __syncwarp();

      int degree = D.degree[candidate];
      int neighOffset = D.offset[candidate];

      for (int j = laneId; j < degree; j += warpSize) {
        int iterBitMask = j / 32; //
        int bitPos = j % 32;
        int neighBitMask =
            levelData.validNeighMask[(start + iter) * maxBitMask + iterBitMask];
        if (neighBitMask & (1U << bitPos)) { //

          ui neigh = D.neighbors[neighOffset + j];

          ui loc = atomicAdd(globalCounter, 1);
          for (int ind = 0; ind < k - 2; ind++) {
            cliqueData.trie[trieSize * ind + loc] =
                levelData.partialCliques[(i) * (k - 1) + ind];
          }
          atomicAdd(&counter[threadIdx.x / warpSize], 1);
          cliqueData.trie[trieSize * (k - 2) + loc] = candidate;
          cliqueData.trie[trieSize * (k - 1) + loc] = neigh;
          cliqueData.status[loc] = -1;
          atomicAdd(&G.cliqueDegree[neigh], 1);
          atomicAdd(&G.cliqueDegree[candidate], 1);
        }
      }
      __syncwarp();

      for (int j = laneId; j < k - 2; j += warpSize) {
        int pClique = levelData.partialCliques[i * (k - 1) + j];
        atomicAdd(&G.cliqueDegree[pClique], counter[threadIdx.x / warpSize]);
      }
      __syncwarp();
    }
  }
}
