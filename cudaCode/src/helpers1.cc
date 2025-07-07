#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cooperative_groups/scan.h>

namespace cg = cooperative_groups;

__global__ void
listIntialCliques(deviceDAGpointer &D, cliqueLevelDataPointer &levelData;
                  ui & label, ui k, ui n, ui m, ui totalWarps, ui totalBlocks) {
  extern __shared__ char sharedMemory[];
  ui sizeOffset = 0;
  ui *counter = (ui *)(sharedMemory + sizeOffset);
  sizeOffset += WARPS_EACH_BLK * sizeof(ui);

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int warpId = idx / warpSize;
  int laneId = idx % warpSize;
  ui partition = warpId;

  cg::thread_block tb = cg::this_thread_block();
  cg::thread_block_tile<warpSize> tile = cg::tiled_partition<warpSize>(tb);

  for (int i = warpId; i < n; i += totalWarps) {
    ui vertex = i;
    ui offset = D.offset[vertex];

    if (laneId == 0) {
      counter[threadId.x / warpSize] = 0;
    }

    __syncwarp();

    ui writeOffset = levelData.offset1[levelData.count[partition]];

    ui degree = D.degree[vertex];
    ui totalActive = 0; // Total active neighbors processed by this warp

    // Process neighbors in chunks of 32
    for (int chunkStart = 0; chunkStart < degree; chunkStart += warpSize) {
      ui chunkSize = min(warpSize, degree - chunkStart);

      // Mask for active threads in this chunk
      ui activeMask = __ballot_sync(0xFFFFFFFF, laneId < chunkSize);

      // Process neighbors in this chunk
      ui loc = 0;
      if (laneId < chunkSize) {
        ui neigh = D.neighbors[offset + chunkStart + laneId];
        ui isActive = (label[neigh] == k);

        // Warp-level prefix sum to calculate the location for each active
        // thread
        loc = cg::exclusive_scan(tile, isActive, cg::plus<ui>());

        if (isActive) {
          label[neigh] = k - 1;
          levelData.candidate1[writeOffset + totalActive + loc] = neigh;
        }
      }

      // Update totalActive for the next chunk
      if (laneId == 0) {
        totalActive += cg::reduce(tile, loc + (laneId < chunkSize ? 1 : 0),
                                  cg::plus<ui>());
      }
    }

    // Update the counter for the warp
    if (laneId == 0) {
      counter[threadId.x / warpSize] += totalActive;
    }

    __syncwarp();

    if (laneId == 0) {
      levelData.offset1[partition] += 1;
    }
  }
}
__global__ void pushRelabel(deviceFlowNetworkPointers flowNetwork,
                            deviceComponentPointers conComp,
                            deviceCliquesPointer finalCliqueData,
                            ui *compCounter, double *upperBound,
                            double *lowerBound, ui *activeNodes,
                            ui *componentsLeft, ui *checkResult,
                            ui *activeCount, ui *syncHelper, ui *shouldBreak,
                            int totalComponents, ui k, ui t, ui partitionSize) {

  int threadId = threadIdx.x;
  int threadsPerBlock = blockDim.x;
  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  for (ui c = 0; c < totalComponents; c++) {
    ui start = conComp.componentOffset[c];
    ui end = conComp.componentOffset[c + 1];
    ui total = end - start;
    double bias = 1.0 / (total * (total - 1));
    bias = max(bias, 0.001);
    if ((upperBound[c] - lowerBound[c]) > bias) {
      ui cliqueStart = compCounter[c];
      ui totalCliques = compCounter[c + 1] - cliqueStart;
      ui fStart = flowNetwork.offset[c];
      ui tFlow = total + totalCliques + 2;
      ui source = tFlow - 2;
      ui sink = tFlow - 1;

      for (ui i = idx; i < tFlow; i += TOTAL_THREAD) {
        flowNetwork.height[fStart + i] = (i == source) ? tFlow : 0;
        flowNetwork.excess[fStart + i] = 0;
        atomicAdd(syncHelper, 1);
      }
      while (true) {
        if (*syncHelper == tFlow) {
          if (idx == 0) {
            *syncHelper = 0;
          }
          break;
        }
      }
      // __syncthreads();

      for (ui j = idx; j < total; j += TOTAL_THREAD) {
        ui nStart = flowNetwork.neighborOffset2[fStart + total + totalCliques];
        ui neighbor = flowNetwork.Edges[nStart + j];
        double capacity = flowNetwork.capacity[nStart + j];

        flowNetwork.flow[nStart + j] = capacity;
        atomicAdd(&flowNetwork.excess[fStart + neighbor], capacity);
        flowNetwork.flow[flowNetwork.neighborOffset2[fStart + neighbor] + 1] =
            -capacity;
        atomicAdd(syncHelper, 1);
      }
      while (true) {
        if (*syncHelper == total) {
          if (idx == 0) {
            *syncHelper = 0;
          }
          break;
        }
      }

      const int maxIterations = 1000;
      for (int iter = 0; iter < maxIterations; iter++) {
        if (idx == 0) {
          *activeCount = 0;
          *syncHelper = 1;
        }
        while (true) {
          if (*syncHelper == 1) {
            if (idx == 0) {
              *syncHelper = 0;
            }
            break;
          }
        }

        // Find active nodes
        for (ui j = idx; j < tFlow; j += TOTAL_THREAD) {
          if (j != source && j != sink &&
              flowNetwork.excess[fStart + j] > 1e-10) {
            ui pos = atomicAdd(&activeCount, 1);
            if (pos < partitionSize) {
              activeNodes[blockIdx.x * partitionSize + pos] = j;
            }
            atomicAdd(syncHelper, 1); // We found work to do
          }
        }
        while (true) {
          if (*syncHelper == tFlow) {
            if (idx == 0) {
              *syncHelper = 0;
            }
            break;
          }
        }

        // Termination check
        if (activeCount == 0)
          break;

        // Process active nodes
        for (ui j = idx; j < activeCount; j += TOTAL_THREAD) {
          ui vertex = activeNodes[blockIdx.x * partitionSize + j];
          ui nStart = flowNetwork.neighborOffset2[fStart + vertex];
          ui nEnd = flowNetwork.neighborOffset2[fStart + vertex + 1];
          bool pushed = false;

          // Try to push flow
          for (ui x = nStart; x < nEnd; x++) {
            ui neighbor = flowNetwork.Edges[x];
            double residual = flowNetwork.capacity[x] - flowNetwork.flow[x];

            if ((flowNetwork.height[fStart + vertex] ==
                 flowNetwork.height[fStart + neighbor] + 1) &&
                (residual > 1e-10)) {

              double delta = min(flowNetwork.excess[fStart + vertex], residual);
              atomicAdd(&flowNetwork.flow[x], delta);
              atomicAdd(&flowNetwork.excess[fStart + vertex], -delta);
              atomicAdd(&flowNetwork.excess[fStart + neighbor], delta);

              // Update backward edge
              ui backStart = flowNetwork.neighborOffset2[fStart + neighbor];
              ui backEnd = flowNetwork.neighborOffset2[fStart + neighbor + 1];
              for (ui y = backStart; y < backEnd; y++) {
                if (flowNetwork.Edges[y] == vertex) {
                  atomicAdd(&flowNetwork.flow[y], -delta);
                  break;
                }
              }
              pushed = true;
            }
          }

          // Relabel if needed
          if (!pushed && flowNetwork.excess[fStart + vertex] > 1e-10) {
            ui minHeight = UINT_MAX;
            for (ui x = nStart; x < nEnd; x++) {
              ui neighbor = flowNetwork.Edges[x];
              double residual = flowNetwork.capacity[x] - flowNetwork.flow[x];
              if (residual > 1e-10) {
                minHeight =
                    min(minHeight, flowNetwork.height[fStart + neighbor]);
              }
            }
            if (minHeight != UINT_MAX) {
              flowNetwork.height[fStart + vertex] = minHeight + 1;
            }
          }
          atomicAdd(syncHelper, 1); // We found work to do
        }

        while (true) {
          if (*syncHelper == *activeCount) {
            if (idx == 0) {
              *syncHelper = 0;
            }
            break;
          }
        }
      }

      if (idx == 0) {
        double alpha = (upperBound[c] + lowerBound[c]) / 2;
        double expectedFlow = totalCliques * k;

        if (fabs(flowNetwork.excess[fStart + sink] - expectedFlow) < 1e-3) {
          upperBound[c] = alpha;
        } else {
          lowerBound[c] = alpha;
          checkResult[c] = 1;
        }

        if ((upperBound[c] - lowerBound[c]) > bias) {
          atomicAdd(componentsLeft, 1);
        }
      }
    }
  }
}