# CDS
## Overview
This repository implements a **GPU-accelerated framework** for the *Densest Subgraph Discovery (DSD)* problem based on **$k$-clique density**.  
It integrates **warp-level $k$-clique enumeration**, **clique-core–based pruning**, and a **parallel push–relabel max-flow solver** to efficiently find subgraphs with the highest higher-order connectivity.

---

## Problem Statement
Given an undirected graph \( G = (V, E) \) and an integer \( k \ge 2 \), the goal is to find a vertex subset \( D \subseteq V \) whose induced subgraph \( G(D) \) maximizes the **$k$-clique density**:

\[
\rho(G(D), \psi) = \frac{n(G(D), \psi)}{|D|}
\]

where \( n(G(D), \psi) \) is the number of $k$-clique instances in \( G(D) \), and each instance represents a fully connected vertex subset of size \( k \).  
This generalizes traditional densest subgraph discovery (based on edges) to **higher-order structures**, enabling the detection of tightly connected communities that go beyond pairwise interactions.

---

##  Build Instructions

**Step 1:** Navigate to the CUDA source directory
```bash
cd cudaCode
```
**Step 2:** Create a build folder and configure
```bash
mkdir -p build
cd build
cmake ..
```

**Step 3:** Complie

```bash
make -j
```

**Tip:** To rebuild cleanly, remove any existing build artifacts first
```bash
rm -rf CMakeFiles CMakeCache.txt
cmake ..
make -j
```
## Run Instructions
You can run the executable using 

```bash
./main <graph_file_path> <k> <pSize> <cpSize> <glBufferSize> <partitionSize> <earlyStop>
```


###  Parameter Description 

| Parameter | Description | Suggested Value | Notes |
|------------|-------------|-----------------|-------|
| `<graph_file_path>` | **Path to the input graph file** (adjacency list format). | -- | Must include the file path in single quotes if it contains spaces. |
| `<k>` | **Clique size** — the order of connectivity to consider. | `2`, `3`, or `4` | Use `2` for edge density, `3` for triangles, `4` for 4-cliques (more complex). |
| `<pSize>` | **Virtual partition size** for storing **partial cliques** during enumeration. | `5000–10000` | Higher values handle larger intermediate sets; increase for denser graphs. |
| `<cpSize>` | **Virtual partition size** for **candidate cliques**. | `10000–20000` | Should typically be ≥ `pSize`; large values reduce re-allocations. |
| `<glBufferSize>` | **Buffer size** for vertices to be removed during **clique-core peeling**. | `50000–100000` | Increase for graphs with large degrees or deep core hierarchies. |
| `<partitionSize>` | **Virtual partition size** for **active nodes** in each **flow network**. | `20000–40000` | Impacts memory layout and parallelism of the push–relabel stage. |
| `<earlyStop>` | **Early stop threshold**. Use `-1` to **disable**. | `-1` | Stops the push relabel algorithm, if the flow doesn't change much in cosecative push relabel loops. This is heuristic ti get approx solution fast |


### Example 
```bash
./main '../../datasets/ca_hepth.txt' 3 8000 16000 60000 30000 -1
```


## Data Format

Your input graph file should be written in **adjacency list format**, using **space-separated integers**.

### File Structure

1. **First line:** 
```bash
<num_vertices> <num_edges>
```
- `num_vertices`: total number of vertices (0-indexed)  
- `num_edges`: total number of **undirected** edges (counted once)

2. **Following lines:**  
Each line represents a vertex and all its neighbors:  
```bash
<vertex_id> <neighbor_1> <neighbor_2> ... <neighbor_t>
```



