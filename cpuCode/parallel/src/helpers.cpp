#include "../inc/helpers.h"
#include "../inc/graph.h"
#include <map>
#include <numeric>
#include <sstream>

CDS ::CDS() {}

CDS ::CDS(Graph *graph, Motif *motif, bool ub1, bool ub2) {
  this->graph = graph;
  this->motif = motif;
  this->ub1 = ub1;
  this->ub2 = ub2;
}

void CDS ::cliqueCoreDecompose(vector<vector<double>> &results) {

  vector<int> mark;
  mark.resize(graph->n, 0);
  vector<int> arrayIndex;
  arrayIndex.resize(graph->n, 0);
  vector<int> newArray;
  newArray.resize(graph->n, 0);
  vector<int> newMap;
  newMap.resize(graph->n, 0);
  // cout << "before clique enum" << endl;

  cliqueEnumerationFast();

  // cout << "after clique enum" << endl;
  unordered_map<int, long> twoDNeighborhood;
  graph->cliqueCore.resize(graph->n, 0);

#pragma omp parallel for
  for (ui i = 0; i < graph->n; i++) {
    graph->cliqueCore[i] = graph->cliqueDegree[i];
  }
  // cout << "here" << endl;

  int totalCliques = 0;

  graph->maxCliqueDegree = 0;

#pragma omp parallel
  {
    ui localMax = 0;
    int localTotal = 0;

#pragma omp for nowait
    for (ui i = 0; i < graph->n; ++i) {
      localTotal += graph->cliqueDegree[i];
      if (graph->cliqueDegree[i] > localMax) {
        localMax = graph->cliqueDegree[i];
      }
    }

#pragma omp critical
    {
      if (localMax > graph->maxCliqueDegree) {
        graph->maxCliqueDegree = localMax;
      }
      totalCliques += localTotal;
    }
  }

  totalCliques = totalCliques / static_cast<double>(motif->size);

  // data structure used to save clique core decompose results

  // cout << "total clique " << totalCliques << endl;

  results.resize(graph->n + 1, vector<double>(5, 0.0));
  results[0][2] = totalCliques / static_cast<double>(graph->n);
  results[0][3] = totalCliques;
  // Create bins for counting the number of vertices with each clique degree
  vector<long> bins;
  bins.resize(graph->maxCliqueDegree + 1, 0);

#pragma omp parallel
  {
    vector<long> localBins(graph->maxCliqueDegree + 1, 0);

#pragma omp for nowait
    for (ui i = 0; i < graph->n; i++) {
      localBins[graph->cliqueDegree[i]]++;
    }

#pragma omp critical
    for (ui i = 0; i <= graph->maxCliqueDegree; i++) {
      bins[i] += localBins[i];
    }
  }

  // cumulative prefix sum of bins
  long sum = 0;
  for (ui i = 0; i <= graph->maxCliqueDegree; i++) {
    long temp = bins[i];
    bins[i] = sum;
    sum += temp;
  }

  // sorts vertices by their clique degree
  vector<ui> pos;
  pos.resize(graph->n, 0);
  vector<ui> sortedVertices;
  sortedVertices.resize(graph->n, 0);

#pragma omp parallel
  {
#pragma omp for
    for (ui i = 0; i < graph->n; i++) {
      ui degree = graph->cliqueDegree[i];
      long position;
#pragma omp atomic capture
      position = bins[degree]++;
      pos[i] = position;
      sortedVertices[position] = i;
    }
  }

  // reset bins to before sort state by shifting one element to right
  for (ui i = graph->maxCliqueDegree; i > 0; i--) {
    bins[i] = bins[i - 1];
  }
  bins[0] = 0;

  int count = 0;

  for (ui i = 0; i < graph->n; ++i) {
    int index = 0;
    long indexMin = 0xFFFFFF;
    for (ui j = 0; j < graph->n; ++j) {
      if (indexMin > graph->cliqueCore[j] && mark[j] == 0) {
        indexMin = graph->cliqueCore[j];
        index = j;
      }
    }
    if (debug) {
      cout << "Removed index: " << index
           << " clique degree: " << graph->cliqueDegree[index] << endl;
    }
    count++;
    results[count][0] = index;
    results[count][1] = graph->cliqueCore[index];
    if (graph->cliqueDegree[index] > 0) {
      // cout << "before get neigh" << endl;

      get2Dneighborhood(twoDNeighborhood, index, mark, arrayIndex, newMap);
      // cout << "after get neigh" << endl;
      long deleteCount = 0;
      if (!twoDNeighborhood.empty()) {
        for (auto &it : twoDNeighborhood) {

          int tempKey = it.first;
          long tempValue = it.second;
          deleteCount += tempValue;
          graph->cliqueCore[tempKey] -= tempValue;
        }
      }
      deleteCount = deleteCount / (motif->size);

      totalCliques -= deleteCount;
      results[count][3] = totalCliques;

      if (graph->n - count > 0) {

        results[count][2] =
            totalCliques / static_cast<double>(graph->n - count);
      } else {
        results[count][2] = 0.0;
      }
    }

    else {
      results[count][3] = totalCliques;
      if (graph->n - count > 0) {

        results[count][2] =
            totalCliques / static_cast<double>(graph->n - count);
      } else {
        results[count][2] = 0.0;
      }
      graph->cliqueCore[index] = 0;
    }
    mark[index] = 1;
  }
}

void CDS::get2Dneighborhood(unordered_map<int, long> &subgraphResults,
                            int index, vector<int> &mark, vector<int> &array,
                            vector<int> &map_s) {
  subgraphResults.clear();
  vector<int> tempList;
  tempList.push_back(index);
  array[index] = 1;
  queue<int> q;
  q.push(index);
  int d = 1;
  while (!q.empty() && (d + 1 <= 2)) {
    int current = q.front();
    q.pop();
    d = array[current];
    for (ui i = 0; i < graph->adjacencyList[current].size(); i++) {
      int neighbor = graph->adjacencyList[current][i];
      if (mark[neighbor] == 0 && array[neighbor] == 0 && (d + 1 <= 2)) {
        array[neighbor] = d + 1;
        tempList.push_back(neighbor);
        q.push(neighbor);
      }
    }
  }
  int count = tempList.size();
  /*cout << "temp list: ";
  for (ui v : tempList) {
    cout << v << " ";
  }
  cout << endl;*/

  vector<int> mapArray;
  mapArray.resize(count, 0);
  int num = 0;
  vector<vector<ui>> subGraph;
  subGraph.resize(count);
  for (int i = 0; i < count; i++) {
    int vertex = tempList[i];
    mapArray[i] = vertex;
    map_s[vertex] = num;
    int tempCount = 0;
    for (int j = 0; j < (int)graph->adjacencyList[vertex].size(); j++) {
      int neighbor = graph->adjacencyList[vertex][j];
      if (array[neighbor] != 0 && neighbor != vertex) {
        subGraph[num].push_back(neighbor);
        tempCount++;
      }
    }
    num++;
  }

  for (int i = 0; i < count; i++) {
    int vertex = tempList[i];
    array[vertex] = 0;
  }
  for (int i = 0; i < count; ++i) {
    // cout << "neigh of " << i << " : ";
    for (int j = 0; j < (int)subGraph[i].size(); ++j) {
      subGraph[i][j] = map_s[subGraph[i][j]];
      // cout << subGraph[i][j] << " ";
    }
    // cout << endl;
  }

  vector<ui> subGraphCliqueDegree;
  subGraphCliqueDegree.resize(count, 0);

  // cout << "clique enum subgraph before" << endl;

  cliqueEnumerationSubgraph(subGraph, subGraphCliqueDegree, motif->size, 0);
  // cout << "clique enum subgraph after" << endl;

  for (int i = 0; i < count; i++) {
    if (subGraphCliqueDegree[i] > 0) {
      subgraphResults[mapArray[i]] = subGraphCliqueDegree[i];
    }
  }
}

void CDS::getlistingOrder(vector<ui> &order) {

  // verticies sorted by reverse core values
  vector<ui> reverseCoreSortedVertices(graph->n);

  coreDecompose(graph->adjacencyList, reverseCoreSortedVertices, graph->degree,
                graph->core, true);
  order.resize(graph->n, 0);
  /*cout << "rev core sorted verticies: ";
  for (ui v : reverseCoreSortedVertices) {
    cout << v << " ";
  }
  cout << endl;*/
  for (ui i = 0; i < reverseCoreSortedVertices.size(); i++) {
    order[reverseCoreSortedVertices[i]] = i + 1;
  }
}

void CDS::coreDecompose(const vector<vector<ui>> adjList,
                        vector<ui> &reverseCoreSortedVertices,
                        vector<ui> &degree, vector<ui> &core, bool fullGraph) {
  ui n = adjList.size();
  reverseCoreSortedVertices.resize(n, 0);
  ui maxDegree = 0;
  core.resize(n, 0);

  for (ui i = 0; i < n; i++) {
    if (degree[i] > maxDegree) {
      maxDegree = degree[i];
    }
    core[i] = degree[i];
  }
  if (fullGraph) {
    graph->maxDegree = maxDegree;
  }

  vector<ui> bins;
  bins.resize(maxDegree + 1, 0);

  for (ui i = 0; i < n; i++) {
    bins[degree[i]]++;
  }

  ui start = 0;
  for (ui i = 0; i <= maxDegree; i++) {
    ui temp = bins[i];
    bins[i] = start;
    start += temp;
  }
  vector<ui> pos;
  pos.resize(n, 0);
  vector<ui> sortedVertices;
  sortedVertices.resize(n, 0);
  for (ui i = 0; i < n; i++) {
    pos[i] = bins[degree[i]];
    sortedVertices[pos[i]] = i;
    bins[degree[i]]++;
  }

  for (ui i = maxDegree; i > 0; i--) {
    bins[i] = bins[i - 1];
  }
  bins[0] = 0;

  for (ui i = 0; i < n; i++) {
    int vertex = sortedVertices[i];
    for (ui j = 0; j < adjList[vertex].size(); j++) {
      int neighbor = adjList[vertex][j];
      if (core[neighbor] > core[vertex]) {
        int du = core[neighbor];
        int pu = pos[neighbor];
        int pw = bins[du];
        int w = sortedVertices[pw];
        if (neighbor != w) {
          pos[neighbor] = pw;
          pos[w] = pu;
          sortedVertices[pu] = w;
          sortedVertices[pw] = neighbor;
        }
        bins[du]++;
        core[neighbor]--;
      }
      reverseCoreSortedVertices[n - i - 1] = vertex;
    }
  }
}
void CDS::cliqueEnumerationFast() {
  vector<ui> order;
  getlistingOrder(order);

  vector<vector<ui>> DAG;
  generateDAG(graph->adjacencyList, DAG, order);

  vector<ui> partialClique;
  vector<ui> candidates;
  for (ui i = 0; i < graph->n; i++) {
    candidates.push_back(i);
  }
  vector<ui> label;
  label.resize(graph->n, motif->size);

  vector<ui> validNeighborCount;
  validNeighborCount.resize(graph->n, 0);
  for (ui i = 0; i < validNeighborCount.size(); i++) {
    validNeighborCount[i] = DAG[i].size();
  }
  graph->cliqueDegree.resize(graph->n, 0);
  graph->totalCliques = 0;

  // Call modified version that parallelizes base case
  listCliquesParallelBase(motif->size, partialClique, candidates, label, DAG,
                          validNeighborCount);
}

void CDS::generateDAG(const vector<vector<ui>> adjList, vector<vector<ui>> &DAG,
                      vector<ui> &order) {
  DAG.resize(adjList.size());
#pragma omp parallel for schedule(static)
  for (ui i = 0; i < adjList.size(); i++) {
    // Phase 1: Count
    int count = 0;
    for (ui j = 0; j < adjList[i].size(); j++) {
      if (order[adjList[i][j]] > order[i])
        count++;
    }

    // Allocate
    DAG[i].resize(count);

    // Phase 2: Fill
    int idx = 0;
    for (ui j = 0; j < adjList[i].size(); j++) {
      if (order[adjList[i][j]] > order[i]) {
        DAG[i][idx++] = adjList[i][j];
      }
    }
  }
}
void CDS::listCliquesParallelBase(ui k, vector<ui> &partialClique,
                                  vector<ui> &candidates, vector<ui> &label,
                                  vector<vector<ui>> &DAG,
                                  vector<ui> &validNeighborCount) {

  if (k == 2) {
    // ========== PARALLEL BASE CASE ==========
    // This is where most work happens - parallelize this!

    // Count cliques from this partial clique
    long cliqueCount = 0;

// Process candidates in parallel
#pragma omp parallel reduction(+ : cliqueCount)
    {

      long localCliqueCount = 0;
#pragma omp for nowait
      for (ui i = 0; i < candidates.size(); i++) {
        ui temp = candidates[i];
        ui bound = validNeighborCount[temp];

        for (ui j = 0; j < bound; j++) {
          localCliqueCount++;
          ui neighbor = DAG[temp][j];
#pragma omp atomic
          graph->cliqueDegree[neighbor]++;
#pragma omp atomic
          graph->cliqueDegree[temp]++;
        }
      }
      cliqueCount += localCliqueCount;
    }
    for (ui v : partialClique) {
      graph->cliqueDegree[v] += cliqueCount;
    }

    graph->totalCliques += cliqueCount;

  } else {
    // ========== RECURSIVE PART - KEEP SERIAL ==========
    // Use your EXACT original recursive code
    for (int i = 0; i < (int)candidates.size(); i++) {
      int temp = candidates[i];
      vector<ui> validNeighbors;

      for (int j = 0; j < (int)DAG[temp].size(); j++) {
        if (label[DAG[temp][j]] == k) {
          label[DAG[temp][j]] = k - 1;
          validNeighbors.push_back(DAG[temp][j]);
        }
      }

      for (int j = 0; j < (int)validNeighbors.size(); j++) {
        ui canTemp = validNeighbors[j];
        int index = 0;

        for (int m = DAG[canTemp].size() - 1; m > index; --m) {
          if (label[DAG[canTemp][m]] == k - 1) {
            while (index < m && label[DAG[canTemp][index]] == k - 1) {
              index++;
            }
            if (label[DAG[canTemp][index]] != k - 1) {
              int temp1 = DAG[canTemp][m];
              DAG[canTemp][m] = DAG[canTemp][index];
              DAG[canTemp][index] = temp1;
            }
          }
        }

        if (DAG[canTemp].size() != 0)
          if (label[DAG[canTemp][index]] == k - 1)
            index++;

        validNeighborCount[canTemp] = index;
      }

      partialClique.push_back(temp);
      listCliquesParallelBase(k - 1, partialClique, validNeighbors, label, DAG,
                              validNeighborCount);
      partialClique.pop_back();

      for (ui j = 0; j < validNeighbors.size(); j++) {
        label[validNeighbors[j]] = k;
      }
    }
  }
}

void CDS::cliqueEnumerationSubgraph(vector<vector<ui>> &subGraph,
                                    vector<ui> &subGraphCliqueDegree,
                                    ui motifSize, ui vertex) {
  vector<ui> reverseCoreSortedVertices(subGraph.size());
  vector<ui> order;
  reverseCoreSortedVertices.resize(subGraph.size(), 0);
  vector<ui> degree;
  degree.resize(subGraph.size(), 0);
  for (ui i = 0; i < subGraph.size(); i++) {
    degree[i] = subGraph[i].size();
  }
  vector<ui> core;
  coreDecompose(subGraph, reverseCoreSortedVertices, degree, core, false);
  // cout << "after core decompose subgraph" << endl;
  order.resize(subGraph.size(), 0);
  for (ui i = 0; i < reverseCoreSortedVertices.size(); i++) {
    order[reverseCoreSortedVertices[i]] = i + 1;
  }

  // for()
  /*cout << "order: ";
  for (ui o : order) {
    cout << o << " ";
  }
  cout << endl;*/
  vector<vector<ui>> DAG;
  // cout << "before DAG: " << endl;
  /*for (ui i = 0; i < subGraph.size(); i++) {
    // cout << "neigh of: " << i << " : ";
    for (ui j = 0; j < subGraph[i].size(); j++) {
      cout << subGraph[i][j] << " ";
    }
    cout << endl;
  }*/
  generateDAG(subGraph, DAG, order);
  /*cout << "after DAG subgraph" << endl;
  cout << DAG.size() << " size of DAG" << endl;
  for (vector<ui> v : DAG) {
    cout << "size " << v.size() << endl;
  }
  for (ui i = 0; i < DAG.size(); i++) {
    cout << "neigh of: " << i << " : ";
    for (ui j = 0; j < DAG[i].size(); j++) {
      cout << DAG[i][j] << " ";
    }
    cout << endl;
  }*/
  vector<ui> partialClique;
  vector<ui> candidates;
  for (ui i = 0; i < subGraph.size(); i++) {
    candidates.push_back(i);
  }
  vector<ui> label;
  label.resize(subGraph.size(), motif->size);
  vector<ui> validNeighborCount;
  validNeighborCount.resize(subGraph.size(), 0);
  for (ui i = 0; i < validNeighborCount.size(); i++) {
    validNeighborCount[i] = DAG[i].size();
  }
  // cout << "before clique list subgraph clique degree " << endl;
  /*for (ui v : subGraphCliqueDegree) {
    cout << v << " ";
  }
  cout << endl;*/

  listCliqueContainsVertex(motifSize, partialClique, candidates, label, DAG,
                           validNeighborCount, subGraphCliqueDegree, vertex);
}

void CDS::listCliqueContainsVertex(ui k, vector<ui> &partialClique,
                                   vector<ui> &candidates, vector<ui> &label,
                                   vector<vector<ui>> &DAG,
                                   vector<ui> &validNeighborCount,

                                   vector<ui> &cliqueDegree, ui vertex) {

  if (k == 2) {
    // cout << "inside k==2" << endl;
    if (debug)
      cout << "----------------------" << endl;
    bool onenode = false;
    string cliqueString = "";
    for (ui i = 0; i < partialClique.size(); i++) {
      cliqueString += to_string(partialClique[i]) + " ";
      if (partialClique[i] == vertex) {
        onenode = true;
      }
    }

    int cliqueCount = 0;
    for (ui i = 0; i < candidates.size(); i++) {

      ui temp = candidates[i];
      for (ui j = 0; j < validNeighborCount[temp]; j++) {
        if (onenode || temp == vertex || DAG[temp][j] == vertex) {
          string wajid = cliqueString + to_string(candidates[i]) + " " +
                         to_string(DAG[temp][j]);
          if (debug)
            cout << wajid << endl;
          cliqueCount++;
          cliqueDegree[DAG[temp][j]]++;
          cliqueDegree[temp]++;
        }
      }
    }

    for (ui i = 0; i < partialClique.size(); i++) {
      int temp = partialClique[i];
      cliqueDegree[temp] += cliqueCount;
    }
    if (debug)
      cout << "----------------------" << endl;

  } else {
    for (int i = 0; i < (int)candidates.size(); i++) {
      int temp = candidates[i];
      vector<ui> validNeighbors;
      // cout << " intial start " << i << " can " << temp << endl;
      for (int j = 0; j < (int)DAG[temp].size(); j++) {
        if (label[DAG[temp][j]] == k) {
          label[DAG[temp][j]] = k - 1;
          validNeighbors.push_back(DAG[temp][j]);
        }
      }
      // cout << " valid neighs size of " << temp << " is "
      //    << validNeighbors.size() << endl;
      for (int j = 0; j < (int)validNeighbors.size(); j++) {

        ui canTemp = validNeighbors[j];
        int index = 0;
        for (int m = DAG[canTemp].size() - 1; m > index; --m) {
          if (label[DAG[canTemp][m]] == k - 1) {
            while (index < m && label[DAG[canTemp][index]] == k - 1) {
              index++;
            }
            if (label[DAG[canTemp][index]] != k - 1) {
              int temp1 = DAG[canTemp][m];
              DAG[canTemp][m] = DAG[canTemp][index];
              DAG[canTemp][index] = temp1;
            }
          }
        }

        if (DAG[canTemp].size() != 0)
          if (label[DAG[canTemp][index]] == k - 1)
            index++;

        validNeighborCount[canTemp] = index;
      }
      partialClique.push_back(temp);
      listCliqueContainsVertex(k - 1, partialClique, validNeighbors, label, DAG,
                               validNeighborCount, cliqueDegree, 0);
      // cout << " pc " << i << " *********" << endl;
      partialClique.pop_back();
      // cout << "here" << endl;
      // cout << validNeighbors.size() << " size " << endl;
      for (ui j = 0; j < validNeighbors.size(); j++) {
        label[validNeighbors[j]] = k;
        // cout << "label " << validNeighbors[j] << endl;
      }
      // cout << "end" << endl;
    }
  }
}

void CDS::locateDensestCore(vector<vector<double>> &coreResults,
                            DensestCoreData &densestCore) {
  graph->maxCliquecore = 0;
  double max = coreResults[0][2];
  for (ui i = 1; i < graph->n; i++) {
    if (max < coreResults[i][2]) {
      // cout << "i: " << coreResults[i][2] << " max: " << max << endl;
      max = coreResults[i][2];
    }
    if (graph->maxCliquecore < coreResults[i][1]) {
      graph->maxCliquecore = coreResults[i][1];
    }
  }

  int lowerBound = (int)ceil(max);

  // int lowerBound = (int)ceil(max);

  int index = 1;
  vector<int> deletedVertices;
  deletedVertices.resize(graph->n, 0);
  for (; index < (int)graph->n; index++) {
    if (coreResults[index][1] >= lowerBound) {
      // deletedVertices.push_back(coreResults[index][0]);
      break;
    }
    deletedVertices[(int)coreResults[index][0]] = -1;
  }

  int temp = 0;
  for (ui i = 0; i < graph->n; i++) {
    if (deletedVertices[i] == 0) {
      deletedVertices[i] = temp;
      temp++;
    }
  }

  vector<vector<double>> newCoreResults;
  newCoreResults.resize(temp, vector<double>(2));

  densestCore.graph.resize(temp);

  ui newGraphSize = temp;

  temp = 0;
  for (ui i = index; i < graph->n; i++) {
    int m = (int)coreResults[i][0];
    newCoreResults[temp][0] = deletedVertices[m];

    newCoreResults[temp][1] = coreResults[i][1];
    temp++;
  }

  for (ui i = 0; i < graph->n; i++) {
    if (deletedVertices[i] != -1) {
      for (ui j = 0; j < graph->adjacencyList[i].size(); j++) {
        int neighbor = graph->adjacencyList[i][j];
        if (deletedVertices[neighbor] != -1) {
          densestCore.graph[deletedVertices[i]].push_back(
              deletedVertices[neighbor]);
        }
      }
    }
  }

  densestCore.reverseMap.resize(newGraphSize, 0);

  for (ui i = 0; i < deletedVertices.size(); i++) {
    if (deletedVertices[i] != -1) {
      densestCore.reverseMap[deletedVertices[i]] = i;
    }
  }
  densestCore.lowerBound = lowerBound;
  densestCore.delVertexIndex = index - 1;
  densestCore.delCliqueCount =
      (int)(coreResults[0][3] - coreResults[index - 1][3]);
  densestCore.density = coreResults[index - 1][2];
  densestCore.maxCliqueCore = graph->maxCliquecore;
}

void CDS::cliqueEnumerationListRecord(
    vector<vector<ui>> newGraph, unordered_map<string, vector<int>> &cliqueData,
    vector<long> &cliqueDegree, ui motifSize) {
  vector<ui> reverseCoreSortedVertices(newGraph.size());
  reverseCoreSortedVertices.resize(newGraph.size(), 0);

  vector<ui> degree;
  degree.resize(newGraph.size(), 0);

#pragma omp parallel for
  for (ui i = 0; i < newGraph.size(); i++) {
    degree[i] = newGraph[i].size();
  }
  vector<ui> core;
  coreDecompose(newGraph, reverseCoreSortedVertices, degree, core, false);

  vector<ui> order;
  order.resize(newGraph.size(), 0);

#pragma omp parallel for
  for (ui i = 0; i < reverseCoreSortedVertices.size(); i++) {
    order[reverseCoreSortedVertices[i]] = i + 1;
  }
  // cout << "after core decom" << endl;
  vector<vector<ui>> DAG;
  generateDAG(newGraph, DAG, order);
  // cout << "after DAG" << endl;

  vector<ui> partialClique;
  vector<ui> candidates;
  for (ui i = 0; i < newGraph.size(); i++) {
    candidates.push_back(i);
  }
  vector<ui> label;
  label.resize(newGraph.size(), motif->size);
  vector<ui> validNeighborCount;
  validNeighborCount.resize(newGraph.size(), 0);
  cliqueDegree.resize(newGraph.size(), 0);

#pragma omp parallel for
  for (ui i = 0; i < validNeighborCount.size(); i++) {
    validNeighborCount[i] = DAG[i].size();
  }

  listCliqueRecordParallelBase(motifSize, partialClique, candidates, label, DAG,
                               validNeighborCount, cliqueData, cliqueDegree);
}

void CDS::listCliqueRecordParallelBase(
    ui k, vector<ui> &partialClique, vector<ui> &candidates, vector<ui> &label,
    vector<vector<ui>> &DAG, vector<ui> &validNeighborCount,
    unordered_map<string, vector<int>> &cliqueData, vector<long> cliqueDegree) {

  if (k == 2) {
    if (debug)
      cout << "----------------------" << endl;
    string cliqueString = "";
    for (ui i = 0; i < partialClique.size(); i++) {
      cliqueString += to_string(partialClique[i]) + " ";
    }

    long cliqueCount = 0;
#pragma omp parallel reduction(+ : cliqueCount)
    {
      long localCliqueCount = 0;

#pragma omp for nowait
      for (ui i = 0; i < candidates.size(); i++) {
        ui temp = candidates[i];
        ui bound = validNeighborCount[temp];

        for (ui j = 0; j < bound; j++) {
          localCliqueCount++;
          ui neighbor = DAG[temp][j];

#pragma omp atomic
          cliqueDegree[neighbor]++;
#pragma omp atomic
          cliqueDegree[temp]++;

          // Build clique string and add to map (critical section)
          string cliqueString1 =
              cliqueString + to_string(temp) + " " + to_string(neighbor);

#pragma omp critical
          {
            vector<int> tempArr(motif->size + 1);
            for (ui x = 0; x < partialClique.size(); x++) {
              tempArr[x] = partialClique[x];
            }
            tempArr[motif->size - 2] = temp;
            tempArr[motif->size - 1] = neighbor;
            tempArr[motif->size] = 1;

            cliqueData[cliqueString1] = tempArr;
          }
        }
      }

      cliqueCount += localCliqueCount;
    }
    for (ui i = 0; i < partialClique.size(); i++) {
      int temp = partialClique[i];
      cliqueDegree[temp] += cliqueCount;
    }

    if (debug)
      cout << "----------------------" << endl;

  } else {
    for (int i = 0; i < (int)candidates.size(); i++) {
      int temp = candidates[i];
      // cout << "temp " << temp << endl;
      vector<ui> validNeighbors;
      for (int j = 0; j < (int)DAG[temp].size(); j++) {
        if (label[DAG[temp][j]] == k) {
          label[DAG[temp][j]] = k - 1;
          validNeighbors.push_back(DAG[temp][j]);
        }
      }
      // cout << validNeighbors.size() << " vn size" << endl;

      for (int j = 0; j < (int)validNeighbors.size(); j++) {

        ui canTemp = validNeighbors[j];
        int index = 0;
        for (int m = DAG[canTemp].size() - 1; m > index; --m) {
          if (label[DAG[canTemp][m]] == k - 1) {
            while (index < m && label[DAG[canTemp][index]] == k - 1) {
              index++;
            }
            if (label[DAG[canTemp][index]] != k - 1) {
              int temp1 = DAG[canTemp][m];
              DAG[canTemp][m] = DAG[canTemp][index];
              DAG[canTemp][index] = temp1;
            }
          }
        }

        if (DAG[canTemp].size() != 0)
          if (label[DAG[canTemp][index]] == k - 1)
            index++;

        validNeighborCount[canTemp] = index;
        // cout << "endhere" << endl;
      }
      partialClique.push_back(temp);
      listCliqueRecordParallelBase(k - 1, partialClique, validNeighbors, label,
                                   DAG, validNeighborCount, cliqueData,
                                   cliqueDegree);
      partialClique.pop_back();
      for (ui j = 0; j < validNeighbors.size(); j++) {
        label[validNeighbors[j]] = k;
      }
    }
  }
}

int CDS::pruneInvalidEdges(vector<vector<ui>> &oldGraph,
                           vector<vector<ui>> &newGraph,
                           unordered_map<string, vector<int>> &cliqueData) {
  // int count = 0;
  vector<unordered_map<int, int>> validEdges(oldGraph.size());

  for (const auto &entry : cliqueData) {
    const vector<int> &temp = entry.second;
    for (ui i = 0; i < temp.size() - 1; i++) {
      for (ui j = 0; j < temp.size() - 1; j++) {
        if (i != j) {
          if (validEdges[temp[i]].find(temp[j]) == validEdges[temp[i]].end()) {
            validEdges[temp[i]][temp[j]] = 0;
          }
        }
      }
    }
  }
  // cout << "here " << endl;

  newGraph.resize(oldGraph.size());

  int totalEdges = 0;

  for (ui i = 0; i < newGraph.size(); i++) {
    for (ui j = 0; j < oldGraph[i].size(); j++) {
      if (validEdges[i].find(oldGraph[i][j]) != validEdges[i].end()) {
        newGraph[i].push_back(oldGraph[i][j]);
        totalEdges += 1;
      }
    }
  }
  // cout << "here 2 " << endl;

  return totalEdges / 2;
}

void CDS::BFS(vector<ui> &status, int vertex, int index,
              const vector<vector<ui>> &newGraph) {
  queue<int> q;
  q.push(vertex);
  while (!q.empty()) {
    int node = q.front();
    q.pop();
    status[node] = index;
    for (int neigh : newGraph[node]) {
      if (status[neigh] == 0) {
        q.push(neigh);
      }
    }
  }
}
void CDS::connectedComponentDecompose(
    vector<vector<ui>> &newGraph,
    unordered_map<string, vector<int>> &cliqueData,
    vector<ConnectedComponentData> &conCompList) {

  // ---------- PHASE 1: Simple Sequential BFS ----------
  vector<ui> status(newGraph.size(), 0);
  int index = 0;
  for (ui i = 0; i < newGraph.size(); i++) {
    if (status[i] == 0) {
      index++;
      BFS(status, i, index, newGraph);
    }
  }

  if (index == 1) {
    ConnectedComponentData conComp;
    conComp.totalCliques = 0;
    conComp.size = newGraph.size();
    conComp.cliqueDegree.resize(conComp.size);
    for (const auto &entry : cliqueData) {
      const vector<int> &temp = entry.second;
      for (ui i = 0; i < temp.size() - 1; i++) {
        conComp.cliqueDegree[temp[i]] += temp[temp.size() - 1];
      }
      conComp.totalCliques += temp[temp.size() - 1];
    }

    conComp.graph = newGraph;
    conComp.cliqueData = cliqueData;
    conComp.density = (double)conComp.totalCliques / ((double)conComp.size);
    conComp.reverseMap.resize(conComp.size);
    iota(conComp.reverseMap.begin(), conComp.reverseMap.end(), 0);
    conCompList.push_back(conComp);

  } else {
    vector<unordered_map<string, vector<int>>> cliqueDataList;
    cliqueDataList.resize(index + 1);
    vector<ui> graphSizeList;
    graphSizeList.resize(index + 1, 0);
    vector<ui> oldToNew;
    oldToNew.resize(newGraph.size(), 0);

    for (ui i = 0; i < newGraph.size(); i++) {
      oldToNew[i] = graphSizeList[status[i]];
      graphSizeList[status[i]]++;
    }

    vector<vector<ui>> reverseMapList;
    reverseMapList.resize(index + 1);

    for (int i = 1; i <= index; i++) {
      reverseMapList.resize(graphSizeList[i]);
    }

    vector<ui> tempIndex;
    tempIndex.resize(index + 1);
    for (ui i = 0; i < newGraph.size(); i++) {
      ui compId = status[i];
      ui newVertexId = oldToNew[i];
      reverseMapList[compId][newVertexId] = i;
      tempIndex[compId]++;
    }

    for (auto &entry : cliqueData) {
      int temp = entry.second[0];
      auto &array = entry.second;
      for (ui i = 0; i < array.size() - 1; i++) {
        array[i] = oldToNew[array[i]];
      }
      cliqueDataList[status[temp]][entry.first] = entry.second;
    }

    vector<vector<vector<ui>>> graphList;
    graphList.resize(index + 1);
    for (int i = 1; i < index + 1; i++) {
      graphList[i].resize(graphSizeList[i]);
    }

    for (ui i = 0; i < newGraph.size(); i++) {
      for (int j = 0; j < (int)newGraph[i].size(); j++) {
        graphList[status[i]][oldToNew[i]].push_back(oldToNew[newGraph[i][j]]);
      }
    }

    // Step 6: Process each component (parallel)
    conCompList.resize(index);
#pragma omp parallel for schedule(dynamic, 1)
    for (int comp = 1; comp <= index; comp++) {
      ConnectedComponentData conComp;
      conComp.size = graphSizeList[comp];
      conComp.reverseMap = reverseMapList[comp];
      conComp.graph = graphList[comp];
      conComp.cliqueData = cliqueDataList[comp];

      // Compute clique degree and total cliques
      conComp.totalCliques = 0;
      conComp.cliqueDegree.resize(conComp.size, 0);

      for (const auto &entry : cliqueDataList[comp]) {
        const vector<int> &temp = entry.second;
        long weight = temp[temp.size() - 1];
        conComp.totalCliques += weight;

        for (ui i = 0; i < temp.size() - 1; i++) {
          ui localVertex = temp[i] - (tempIndex[comp] - graphSizeList[comp]);
          conComp.cliqueDegree[localVertex] += weight;
        }
      }

      conComp.density = (conComp.size > 0) ? (double)conComp.totalCliques /
                                                 (double)conComp.size
                                           : 0.0;

      conCompList[comp - 1] = conComp;
    }
  }
}
void CDS::dynamicExact(vector<ConnectedComponentData> &conCompList,
                       DensestCoreData &densestCore,
                       finalResult &densestSubgraph, bool ub1, bool ub2) {

  // Sequential initialization (keep as is)
  ConnectedComponentData current = conCompList[0];
  float lowerBound = 0.0f;

  for (ui i = 0; i < conCompList.size(); i++) {
    if (lowerBound < conCompList[i].density) {
      lowerBound = conCompList[i].density;
      current = conCompList[i];
    }
  }

  if (ceil(lowerBound) < ceil(densestCore.density)) {
    lowerBound = densestCore.density;
  }

  float upperBound = densestCore.maxCliqueCore;

  // Initialize densest subgraph
  densestSubgraph.verticies.resize(current.size);

#pragma omp parallel for
  for (int i = 0; i < current.size; i++) {
    densestSubgraph.verticies[i] =
        densestCore.reverseMap[current.reverseMap[i]];
  }
  densestSubgraph.size = current.size;
  densestSubgraph.density = lowerBound;

  // Iterate over components (sequential)
  vector<int> res;

  for (ui cid = 0; cid < conCompList.size(); cid++) {
    current = conCompList[cid];

    // UB1 computation (small, keep sequential)
    if (ub1) {
      float k = (float)motif->size;
      float log_fact = 0.0f;
      for (ui i = 1; i <= motif->size; ++i)
        log_fact += logf((float)i);

      float dem = expf(log_fact / k);
      float num = powf((float)current.totalCliques, (k - 1.0f) / k);
      float ub1_val = num / dem;

      if (ub1_val < upperBound)
        upperBound = ub1_val;
    }

    // UB2 computation - PARALLELIZE THIS (most expensive part of bounds)
    if (ub2) {
      float dem = (float)(motif->size * current.totalCliques);
      float ub2_val = 0.0f;

// PARALLEL reduction
#pragma omp parallel for reduction(+ : ub2_val)
      for (ui v = 0; v < current.cliqueDegree.size(); v++) {
        float pv = current.cliqueDegree[v] / dem;
        ub2_val += pv * pv;
      }

      ub2_val *= current.totalCliques;
      if (ub2_val < upperBound)
        upperBound = ub2_val;
    }

    // Exact solve (sequential - it's max-flow)
    exact(res, current, lowerBound, upperBound);

    // Count motifs - PARALLELIZE THIS
    long cliqueCount = 0;

    // Convert cliqueData to vector for safe parallel access
    vector<pair<string, vector<int>>> cliqueVec(current.cliqueData.begin(),
                                                current.cliqueData.end());

#pragma omp parallel for reduction(+ : cliqueCount)
    for (ui idx = 0; idx < cliqueVec.size(); idx++) {
      const vector<int> &temp = cliqueVec[idx].second;
      bool allInResult = true;

      for (ui i = 0; i < temp.size() - 1; ++i) {
        if (res[temp[i]] == -1) {
          allInResult = false;
          break;
        }
      }

      if (allInResult) {
        cliqueCount += temp[temp.size() - 1]; // Add the count
      }
    }

    // Count vertices - PARALLELIZE THIS
    ui vertexCount = 0;
#pragma omp parallel for reduction(+ : vertexCount)
    for (int i = 0; i < current.size; i++) {
      if (res[i] != -1)
        vertexCount++;
    }

    if (vertexCount == 0)
      vertexCount = current.size;

    float density = (float)cliqueCount / (float)vertexCount;

    // Update densest subgraph (sequential - need to update vectors)
    if (density > densestSubgraph.density) {
      densestSubgraph.density = density;
      lowerBound = density;
      densestSubgraph.size = vertexCount;

      densestSubgraph.verticies.clear();
      densestSubgraph.verticies.reserve(vertexCount);

      // Build vertices list (could parallelize but it's small)
      for (int i = 0; i < current.size; i++) {
        if (res[i] != -1) {
          densestSubgraph.verticies.push_back(
              densestCore.reverseMap[current.reverseMap[i]]);
        }
      }
    }
  }
}
void CDS::exact(vector<int> &res, ConnectedComponentData &C, float lowerBound,
                float upperBound) {

  float alpha = (lowerBound + upperBound) * 0.5f;
  float bias = 1.0f / (C.size * (C.size - 1));
  if (bias < 1e-2f)
    bias = 1e-2f;

  FlowNetwork FN;
  createFlownetwork(FN, C, alpha);
  /*cout << "create flow network" << endl;

  for (ui i = 0; i < FN.G.size(); i++) {
    cout << "Edges of node " << i << ": ";
    for (auto &e : FN.G[i]) {
      cout << e.to << " ";
    }
    cout << endl;
  }*/

  vector<pair<int, int>> parent(FN.G.size());
  res.assign(C.size, 1);

  float target = C.totalCliques * motif->size;

  while ((upperBound - lowerBound) > bias) {

    float flow = edmondsKarp(FN, parent);
    // cout << "create edmond" << endl;

    if (fabs(flow - target) < 1e-2f) {
      upperBound = alpha;
    } else {
      lowerBound = alpha;

      // Java-equivalent: res[i] = parent[i] != -1 ? parent[i] : -1
      for (int i = 0; i < C.size; i++)
        res[i] = parent[i].first;
    }

    alpha = (lowerBound + upperBound) * 0.5f;
    updateFlownetwork(FN, C, alpha);
    // /cout << "create update network" << endl;
  }
}

void CDS::createFlownetwork(FlowNetwork &FN, ConnectedComponentData &C,
                            float alpha) {

  int totalNodes = C.size + C.totalCliques + 2;
  FN.init(totalNodes);

  FN.source = C.size + C.totalCliques;
  FN.sink = FN.source + 1;

// Reserve memory in parallel
#pragma omp parallel for
  for (int i = 0; i < totalNodes; i++) {
    if (i < C.size)
      FN.G[i].reserve(C.cliqueDegree[i] + 2);
    else if (i == FN.source || i == FN.sink)
      FN.G[i].reserve(C.size);
    else
      FN.G[i].reserve(motif->size * 2);
  }

  // Convert cliqueData to vector for parallel processing
  vector<vector<int>> cliqueVec;
  cliqueVec.reserve(C.cliqueData.size());
  for (const auto &entry : C.cliqueData) {
    cliqueVec.push_back(entry.second);
  }

// Add source->vertex and vertex->sink edges in parallel
#pragma omp parallel for
  for (int v = 0; v < C.size; v++) {
#pragma omp critical
    {
      FN.addEdge(FN.source, v, (float)C.cliqueDegree[v]);
      FN.addEdge(v, FN.sink, alpha * motif->size);
      FN.sinkEdgeIdx[v] = FN.G[v].size() - 1;
    }
  }

// Add clique edges in parallel
#pragma omp parallel for
  for (ui idx = 0; idx < cliqueVec.size(); idx++) {
    const vector<int> &clq = cliqueVec[idx];
    float count = clq[motif->size];
    float w = count * (motif->size - 1);
    int cliqueIdx = C.size + idx;

// Need critical section for addEdge or pre-allocate
#pragma omp critical
    {
      for (ui i = 0; i < motif->size; i++) {
        int v = clq[i];
        FN.addEdge(cliqueIdx, v, w);
        FN.addEdge(v, cliqueIdx, count);
      }
    }
  }
}
float CDS::edmondsKarp(FlowNetwork &FN, vector<pair<int, int>> &parent) {
  float sum = 0.0f;
  float pushed;

  while ((pushed = augmentPath(FN, parent)) != -1) {
    // Augment path updates (sequential)
    for (int v = FN.sink; v != FN.source;) {
      auto [u, ei] = parent[v];
      Edge &e = FN.G[u][ei];
      Edge &rev = FN.G[v][e.rev];

      e.cap -= pushed;
      rev.cap += pushed;

      v = u;
    }
    sum += pushed;
  }
  return sum;
}

float CDS::augmentPath(FlowNetwork &FN, vector<pair<int, int>> &parent) {
  // BFS is inherently sequential, but we can parallelize frontier expansion
  fill(parent.begin(), parent.end(), make_pair(-1, -1));

  static vector<int> q;
  q.clear();
  q.push_back(FN.source);
  parent[FN.source] = {FN.source, -1};

  // Sequential BFS (tried parallel BFS, not worth the complexity)
  for (int qi = 0; qi < (int)q.size(); qi++) {
    int u = q[qi];
    if (u == FN.sink)
      break;

    // Process neighbors - could parallelize but frontier is usually small
    for (int ei = 0; ei < (int)FN.G[u].size(); ei++) {
      Edge &e = FN.G[u][ei];
      if (parent[e.to].first == -1 && e.cap > 0) {
        parent[e.to] = {u, ei};
        q.push_back(e.to);
      }
    }
  }

  if (parent[FN.sink].first == -1)
    return -1;

  // Find bottleneck (sequential)
  float flow = FLT_MAX;
  for (int v = FN.sink; v != FN.source;) {
    auto [u, ei] = parent[v];
    flow = min(flow, FN.G[u][ei].cap);
    v = u;
  }

  return flow;
}
void CDS::updateFlownetwork(FlowNetwork &FN, ConnectedComponentData &C,
                            float alpha) {

// Reset all capacities in parallel
#pragma omp parallel for
  for (size_t i = 0; i < FN.G.size(); i++) {
    auto &u = FN.G[i];
    for (auto &e : u) {
      e.cap = e.max;
    }
  }

  // Update vertex → sink edges in parallel
  const int motif_size = motif->size; // Cache outside parallel region
#pragma omp parallel for
  for (int v = 0; v < C.size; v++) {
    int idx = FN.sinkEdgeIdx[v];
    FN.G[v][idx].cap = FN.G[v][idx].max = alpha * motif_size;
  }
}

void CDS::DSD() {
  auto start = Clock::now();
  vector<vector<double>> results;
  auto t1 = Clock::now();
  cliqueCoreDecompose(results);
  auto t2 = Clock::now();
  double time_cd = std::chrono::duration<double, std::milli>(t2 - t1).count();

  if (debug) {
    for (ui i = 0; i < results.size(); i++) {
      cout << results[i][0] << " | " << results[i][1] << " | " << results[i][2]
           << " |  " << results[i][3] << " | " << results[i][4] << endl;
    }
  }

  DensestCoreData densestCore;
  t1 = Clock::now();

  locateDensestCore(results, densestCore);
  t2 = Clock::now();
  double time_lc = std::chrono::duration<double, std::milli>(t2 - t1).count();
  if (debug) {
    cout << "LowerBound: " << densestCore.lowerBound << endl;
    cout << "Del vertex index: " << densestCore.delVertexIndex << endl;
    cout << "Del clique count: " << densestCore.delCliqueCount << endl;
    cout << "Densest Core Density: " << densestCore.density << endl;
    cout << "Max Clique core value: " << densestCore.maxCliqueCore << endl;
    for (ui i = 0; i < densestCore.graph.size(); i++) {
      cout << "Neigh of: " << densestCore.reverseMap[i] << " : ";
      for (ui j = 0; j < densestCore.graph[i].size(); j++) {
        cout << densestCore.reverseMap[densestCore.graph[i][j]] << " ";
      }
      cout << endl;
    }
  }

  unordered_map<string, vector<int>> cliqueData;
  vector<long> cliqueDegree;
  t1 = Clock::now();

  cliqueEnumerationListRecord(densestCore.graph, cliqueData, cliqueDegree,
                              motif->size);
  t2 = Clock::now();
  double time_cl = std::chrono::duration<double, std::milli>(t2 - t1).count();
  if (debug) {
    cout << "Cliques " << endl;
    for (const auto &entry : cliqueData) {
      cout << "Key : " << entry.first << endl;
      const vector<int> &temp = entry.second;
      for (ui i = 0; i < temp.size() - 1; i++) {
        cout << temp[i] << " ";
      }
      cout << endl;
    }
  }

  vector<vector<ui>> newGraph;
  t1 = Clock::now();

  int validEdgeCount =
      pruneInvalidEdges(densestCore.graph, newGraph, cliqueData);
  t2 = Clock::now();
  double time_pe = std::chrono::duration<double, std::milli>(t2 - t1).count();
  if (debug) {
    cout << "Graph after pruning:  " << endl;
    for (ui i = 0; i < newGraph.size(); i++) {
      cout << "Neigh of i: " << i << " : ";
      for (ui j = 0; j < newGraph[i].size(); j++) {
        cout << newGraph[i][j] << " ";
      }
      cout << endl;
    }
    cout << "New Neighbor size: " << validEdgeCount << endl;
  }

  vector<ConnectedComponentData> conCompList;
  t1 = Clock::now();

  connectedComponentDecompose(newGraph, cliqueData, conCompList);
  t2 = Clock::now();
  double time_cc = std::chrono::duration<double, std::milli>(t2 - t1).count();
  if (debug) {
    cout << "Connected Components : " << endl;

    for (ConnectedComponentData cc : conCompList) {
      cout << "Size: " << cc.size << endl;
      cout << "Clique Count: " << cc.totalCliques << endl;
      cout << "Density: " << cc.density << endl;
      for (ui i = 0; i < cc.graph.size(); i++) {
        cout << "Neigh of i: " << i << " : ";
        for (ui j = 0; j < cc.graph[i].size(); j++) {
          cout << cc.graph[i][j] << " ";
        }
        cout << endl;
      }
      for (const auto &entry : cc.cliqueData) {
        cout << "Key : " << entry.first << endl;
        const vector<int> &temp = entry.second;
        for (ui i = 0; i < temp.size() - 1; i++) {
          cout << temp[i] << " ";
        }
        cout << endl;
      }
    }
  }

  finalResult densestSubgraph;
  t1 = Clock::now();
  float ub1_val;
  float ub2_val;
  dynamicExact(conCompList, densestCore, densestSubgraph, ub1_val, ub2_val, ub1,
               ub2);
  cout << "Dynamic exact: " << endl;
  t2 = Clock::now();
  double time_de = std::chrono::duration<double, std::milli>(t2 - t1).count();
  auto end = Clock::now();

  double time_ms =
      std::chrono::duration<double, std::milli>(end - start).count();

  cout << "Execution_time: " << time_ms << " ms" << endl;
  cout << "Component_decompose: " << time_cd << " ms" << endl;
  cout << "Locate_core: " << time_lc << " ms" << endl;
  cout << "Clique_Listing: " << time_cl << " ms" << endl;
  cout << "Prune_Edges: " << time_pe << " ms" << endl;
  cout << "Connected_components: " << time_cc << " ms" << endl;
  cout << "Dynamic_Exact: " << time_de << " ms" << endl;

  cout << "Total_cliques: " << graph->totalCliques << endl;
  cout << "Densest_core_deleted_cliques: " << densestCore.delCliqueCount
       << endl;
  cout << "Denseset_core_kmax: " << densestCore.maxCliqueCore << endl;
  cout << "Denseset_core_density: " << densestCore.density << endl;
  cout << "Densest_core_size: " << densestCore.graph.size() << endl;
  int edgeCount = 0;
  for (vector<ui> v : densestCore.graph) {
    edgeCount += v.size();
  }
  cout << "Densest_core_edges_count: " << edgeCount / 2 << endl;
  cout << "Remaining_edges_after_prune: " << validEdgeCount << endl;
  cout << "Component_number: " << conCompList.size() << endl;

  double lb = 0;
  for (ui i = 0; i < conCompList.size(); i++) {
    if (lb < conCompList[i].density) {
      lb = conCompList[i].density;
    }
  }
  cout << "Paper_lower_bound: " << lb << endl;
  cout << "Paper_upper_bopund: " << densestCore.maxCliqueCore << endl;

  cout << "K: " << motif->size << endl;
  cout << "Density: " << densestSubgraph.density << endl;
  cout << "Size: " << densestSubgraph.size << endl;
  cout << "densest subgraph: ";

  for (ui v : densestSubgraph.verticies) {
    cout << v << " ";
  }
  cout << endl;
}
