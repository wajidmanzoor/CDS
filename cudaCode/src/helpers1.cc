#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#include <cooperative_groups/scan.h>

namespace cg = cooperative_groups;

__global__ void listIntialCliques(deviceDAGpointer &D, cliqueLevelDataPointer &levelData; ui &label,ui k, ui n, ui m, ui totalWarps, ui totalBlocks) {
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
            counter[threadId.x/warpSize] = 0;
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

                // Warp-level prefix sum to calculate the location for each active thread
                loc = cg::exclusive_scan(tile, isActive, cg::plus<ui>());

                if (isActive) {
                    label[neigh] = k - 1;
                    levelData.candidate1[writeOffset + totalActive + loc] = neigh;
                }
            }

            // Update totalActive for the next chunk
            if (laneId == 0) {
                totalActive += cg::reduce(tile, loc + (laneId < chunkSize ? 1 : 0), cg::plus<ui>());
            }
        }

        // Update the counter for the warp
        if (laneId == 0) {
            counter[threadId.x/warpSize] += totalActive;
        }

        __syncwarp();

        if (laneId == 0) {
            levelData.offset1[partition] += 1;
        }
    }

    
}