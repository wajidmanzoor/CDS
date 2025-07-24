__global__ void listIntialCliques(deviceDAGpointer D,
                                  cliqueLevelDataPointer levelData, ui *label,
                                  ui k, ui n, ui psize, ui cpSize,
                                  ui maxBitMask, ui level, ui totalWarps,
                                  size_t partialSize, size_t candidateSize,
                                  size_t maskSize, size_t offsetSize) {
  extern __shared__ char sharedMemory[];
  ui *counter = (ui *)sharedMemory;

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int warpId = idx / warpSize;
  int laneId = idx % warpSize;

  size_t oP = (size_t)warpId * psize;
  size_t oCP = (size_t)warpId * cpSize;
  size_t oO = (size_t)warpId * ((psize / (k - 1)) + 1);
  size_t oM = oCP * maxBitMask;
  size_t bitCap = ((size_t)TOTAL_WARPS * n + 7) >> 3;

  for (ui i = warpId; i < n; i += totalWarps) {
    ui vertex = i, vOff = D.offset[vertex];

    if (laneId == 0)
      counter[threadIdx.x / warpSize] = 0;
    __syncwarp();

    size_t cOff =
        oCP + levelData.offsetPartition[oO + levelData.count[warpId + 1]];

    for (ui j = laneId; j < D.degree[vertex]; j += warpSize) {
      ui neigh = D.neighbors[vOff + j];
      size_t bit = (size_t)warpId * n + neigh;
      size_t bIdx = bit / 8;
      uint8_t mask = 1 << (bit % 8);

      uint8_t old = atomicOr(&label[bIdx], mask);
      if (!(old & mask)) {
        ui loc = atomicAdd(&counter[threadIdx.x / warpSize], 1);

        levelData.candidatesPartition[cOff + loc] = neigh;
      }
    }

    __syncwarp();

    if (laneId == 0 && counter[threadIdx.x / warpSize] > 0) {
      size_t w = oP + levelData.count[warpId + 1] * (k - 1) + level;

      levelData.partialCliquesPartition[w] = vertex;

      levelData.count[warpId + 1] += 1;

      levelData.offsetPartition[oO + levelData.count[warpId + 1]] =
          levelData.offsetPartition[oO + levelData.count[warpId + 1] - 1] +
          counter[threadIdx.x / warpSize];
    }

    __syncwarp();

    for (ui j = laneId; j < counter[threadIdx.x / warpSize]; j += warpSize) {
      ui cand = levelData.candidatesPartition[cOff + j];
      ui dOff = D.offset[cand], deg = D.degree[cand], chunks = (deg + 31) / 32;

      for (ui m = 0; m < chunks; m++) {
        ui bitmask = 0, from = m * 32, to = min(from + 32, deg);
        for (ui x = from; x < to; x++) {
          ui nb = D.neighbors[dOff + x];
          size_t bit = (size_t)warpId * n + nb;
          size_t bIdx = bit / 8;
          uint8_t msk = 1 << (bit % 8);

          if ((label[bIdx] & msk) != 0)
            bitmask |= 1 << (x - from);
        }

        size_t maskIdx =
            oM +
            (levelData.offsetPartition[oO + levelData.count[warpId + 1] - 1] +
             j) *
                maxBitMask +
            m;

        levelData.validNeighMaskPartition[maskIdx] = bitmask;
      }
    }

    __syncwarp();

    for (ui x = laneId; x < n; x += warpSize) {
      size_t bit = (size_t)warpId * n + x;
      size_t bIdx = bit / 8;
      uint8_t msk = ~(1 << (bit % 8));
      atomicAnd(&label[bIdx], msk);
    }

    __syncwarp();
  }
}
