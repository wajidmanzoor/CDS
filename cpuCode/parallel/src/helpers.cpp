#include "../inc/helpers.h"
#include <map>
#include <numeric>
#include <sstream>

CDS ::CDS() {}

CDS ::CDS(Graph *graph, Motif *motif) {
  this->graph = graph;
  this->motif = motif;
}

void CDS ::cliqueCoreDecompose(vector<vector<double>> &results) {

  vector<int> mark;
  mark.resize(graph->n, 0);
  vector<int> arrayIndex;
  arrayIndex.resize(graph->n, 0);
  vector<int> newArray;
  arrayIndex.resize(graph->n, 0);
  vector<int> newMap;
  newMap.resize(graph->n, 0);

  cliqueEnumerationFast();

  unordered_map<int, long> twoDNeighborhood;
#pragma omp parallel for schedule(static)
  for (ui i = 0; i < graph->n; i++) {
    graph->cliqueCore[i] = graph->cliqueDegree[i];
  }

  int totalCliques = 0;
  int maxCliqueDegreeLocal = 0;

#pragma omp parallel for reduction(+ : totalCliques)                           \
    reduction(max : maxCliqueDegreeLocal)
  for (ui i = 0; i < graph->n; i++) {
    maxCliqueDegreeLocal =
        std::max(maxCliqueDegreeLocal, (int)graph->cliqueDegree[i]);
    totalCliques += graph->cliqueDegree[i];
  }

  // write back once (sequential, safe)
  graph->maxCliqueDegree = maxCliqueDegreeLocal;

  totalCliques = totalCliques / static_cast<double>(motif->size);

  // data structure used to save clique core decompose results

  results.resize(graph->n + 1, vector<double>(5, 0.0));
  results[0][2] = totalCliques / static_cast<double>(graph->n);
  results[0][3] = totalCliques;
  // Create bins for counting the number of vertices with each clique degree
  vector<long> bins;
  bins.resize(graph->maxCliqueDegree + 1, 0);

  int maxD = graph->maxCliqueDegree;
  int T = omp_get_max_threads();

  vector<vector<long>> localBins(T, vector<long>(maxD + 1, 0));

#pragma omp parallel
  {
    int tid = omp_get_thread_num();
#pragma omp for
    for (ui i = 0; i < graph->n; i++) {
      localBins[tid][graph->cliqueDegree[i]]++;
    }
  }

  for (int t = 0; t < T; t++) {
    for (int d = 0; d <= maxD; d++) {
      bins[d] += localBins[t][d];
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
  for (ui i = 0; i < graph->n; i++) {
    pos[i] = bins[graph->cliqueDegree[i]];
    sortedVertices[pos[i]] = i;
    bins[graph->cliqueDegree[i]]++;
  }

  // reset bins to before sort state by shifting one element to right
  for (ui i = graph->maxCliqueDegree; i > 0; i--) {
    bins[i] = bins[i - 1];
  }
  bins[0] = 0;

  int count = 0;

  for (ui i = 0; i < graph->n; i++) {
    int index = -1;
    long indexMin = LONG_MAX;

#pragma omp parallel
    {
      int localIndex = -1;
      long localMin = LONG_MAX;

#pragma omp for nowait
      for (ui j = 0; j < graph->n; j++) {
        if (!mark[j] && graph->cliqueDegree[j] < localMin) {
          localMin = graph->cliqueDegree[j];
          localIndex = j;
        }
      }

#pragma omp critical
      {
        if (localMin < indexMin) {
          indexMin = localMin;
          index = localIndex;
        }
      }
    }
    if (debug) {
      cout << "Removed index: " << index
           << " clique degree: " << graph->cliqueDegree[index] << endl;
    }
    count++;
    results[count][0] = index;
    results[count][1] = graph->cliqueDegree[index];
    if (graph->cliqueDegree[index] > 0) {
      get2Dneighborhood(twoDNeighborhood, index, mark, arrayIndex, newMap);
      long deleteCount = 0;
      if (!twoDNeighborhood.empty()) {
        for (auto &it : twoDNeighborhood) {

          int tempKey = it.first;
          long tempValue = it.second;
          deleteCount += tempValue;
          graph->cliqueCore[tempKey] -= tempValue;
        }
        deleteCount = deleteCount / static_cast<double>(motif->size - 1);
        totalCliques -= deleteCount;
        results[count][3] = totalCliques;
        if (graph->n - count > 0) {
          results[count][2] =
              totalCliques / static_cast<double>(graph->n - count);
        } else {
          results[count][2] = 0.0;
        }
      }
    } else {
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

  vector<int> mapArray;
  mapArray.resize(count, 0);
  int num = 0;
  vector<vector<ui>> subGraph;
  subGraph.resize(count);
  for (ui i = 0; i < count; i++) {
    int vertex = tempList[i];
    mapArray[i] = vertex;
    map_s[vertex] = num;
    int tempCount = 0;
    for (int j = 0; j < graph->adjacencyList[vertex].size(); j++) {
      int neighbor = graph->adjacencyList[vertex][j];
      if (array[neighbor] == 0 && neighbor != vertex) {
        subGraph[num].push_back(neighbor);
        tempCount++;
      }
    }
    num++;
  }

  for (ui i = 0; i < count; i++) {
    int vertex = tempList[i];
    array[vertex] = 0;
  }
  for (int i = 0; i < count; ++i) {
    for (int j = 0; j < subGraph[i].size(); ++j) {
      subGraph[i][j] = map_s[subGraph[i][j]];
    }
  }

  vector<ui> subGraphCliqueDegree;
  subGraphCliqueDegree.resize(count, 0);

  cliqueEnumerationSubgraph(subGraph, subGraphCliqueDegree, motif->size, 0);

  for (ui i = 0; i < count; i++) {
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
    core[vertex] = degree[vertex];
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

  listCliques(motif->size, partialClique, candidates, label, DAG,
              validNeighborCount);
}

void CDS::generateDAG(const vector<vector<ui>> adjList, vector<vector<ui>> &DAG,
                      vector<ui> &order) {
  int count;
  DAG.resize(adjList.size());
#pragma omp parallel for schedule(dynamic)
  for (ui i = 0; i < adjList.size(); i++) {
    int count = 0;
    for (ui j = 0; j < adjList[i].size(); j++) {
      if (order[adjList[i][j]] > order[i])
        count++;
    }
    DAG[i].resize(count);

    int idx = 0;
    for (ui j = 0; j < adjList[i].size(); j++) {
      if (order[adjList[i][j]] > order[i])
        DAG[i][idx++] = adjList[i][j];
    }
  }
}

void CDS::listCliques(ui k, vector<ui> &partialClique, vector<ui> &candidates,
                      vector<ui> &label, vector<vector<ui>> &DAG,
                      vector<ui> &validNeighborCount) {

  if (k == 2) {
    string cliqueString = "";
    for (ui i = 0; i < partialClique.size(); i++) {
      cliqueString += to_string(partialClique[0]) + " ";
    }

    int cliqueCount = 0;
    for (ui i = 0; i < candidates.size(); i++) {
      int temp = candidates[i];
      for (int j = 0; j < validNeighborCount[temp]; j++) {
        cliqueCount++;
        graph->totalCliques++;
        graph->cliqueDegree[DAG[temp][j]]++;
        graph->cliqueDegree[temp]++;
      }
    }

    for (ui i = 0; i < candidates.size(); i++) {
      int temp = candidates[i];
      graph->cliqueDegree[temp] += cliqueCount;
    }

  } else {
    for (int i = 0; i < candidates.size(); i++) {
      int temp = candidates[i];
      vector<ui> validNeighbors;
      for (int j = 0; j < DAG[temp].size(); j++) {
        if (label[DAG[temp][j]] == k) {
          label[DAG[temp][j]] = k - 1;
          validNeighbors.push_back(DAG[temp][j]);
        }
      }

      for (int j = 0; j < validNeighbors.size(); j++) {

        ui canTemp = validNeighbors[j];
        int index = 0;
        for (ui m = DAG[canTemp].size() - 1; m > index; --m) {
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

        if (DAG[canTemp].size() != 0 && label[DAG[canTemp][index]] == k - 1)
          index++;

        validNeighborCount[canTemp] = index;
      }
      partialClique.push_back(temp);
      listCliques(k - 1, partialClique, validNeighbors, label, DAG,
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

  order.resize(subGraph.size(), 0);
  for (ui i = 0; i < reverseCoreSortedVertices.size(); i++) {
    order[reverseCoreSortedVertices[i]] = i + 1;
  }
  vector<vector<ui>> DAG;
  generateDAG(subGraph, DAG, order);
  vector<ui> partialClique;
  vector<ui> candidates;
  for (ui i = 0; i < subGraph.size(); i++) {
    candidates.push_back(i);
  }
  vector<ui> label;
  label.resize(subGraph.size(), motif->size);
  vector<ui> validNeighborCount;
  validNeighborCount.resize(subGraph.size(), 0);

  listCliqueContainsVertex(motif->size, partialClique, candidates, label, DAG,
                           validNeighborCount, subGraphCliqueDegree, 0);
}

void CDS::listCliqueContainsVertex(ui k, vector<ui> &partialClique,
                                   vector<ui> &candidates, vector<ui> &label,
                                   vector<vector<ui>> &DAG,
                                   vector<ui> &validNeighborCount,
                                   vector<ui> &cliqueDegree, ui vertex) {
  if (k == 2) {
    bool onenode = false;
    string cliqueString = "";
    for (ui i = 0; i < partialClique.size(); i++) {
      cliqueString += to_string(partialClique[0]) + " ";
      if (partialClique[i] == vertex) {
        onenode = true;
      }
    }

    int cliqueCount = 0;
    for (ui i = 0; i < candidates.size(); i++) {
      int temp = candidates[i];
      for (int j = 0; j < validNeighborCount[temp]; j++) {
        if (onenode || temp == vertex || DAG[temp][j] == vertex) {
          cliqueCount++;
          cliqueDegree[DAG[temp][j]]++;
          cliqueDegree[temp]++;
        }
      }
    }

    for (ui i = 0; i < candidates.size(); i++) {
      int temp = candidates[i];
      cliqueDegree[temp] += cliqueCount;
    }

  } else {
    for (int i = 0; i < candidates.size(); i++) {
      int temp = candidates[i];
      vector<ui> validNeighbors;
      for (int j = 0; j < DAG[temp].size(); j++) {
        if (label[DAG[temp][j]] == k) {
          label[DAG[temp][j]] = k - 1;
          validNeighbors.push_back(DAG[temp][j]);
        }
      }

      for (int j = 0; j < validNeighbors.size(); j++) {

        ui canTemp = validNeighbors[j];
        int index = 0;
        for (ui m = DAG[canTemp].size() - 1; m > index; --m) {
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

        if (DAG[canTemp].size() != 0 && label[DAG[canTemp][index]] == k - 1)
          index++;

        validNeighborCount[canTemp] = index;
      }
      partialClique.push_back(temp);
      listCliqueContainsVertex(motif->size, partialClique, candidates, label,
                               DAG, validNeighborCount, cliqueDegree, 0);
      partialClique.pop_back();
      for (ui j = 0; j < validNeighbors.size(); j++) {
        label[validNeighbors[j]] = k;
      }
    }
  }
}

void CDS::locateDensestCore(vector<vector<double>> &coreResults,
                            DensestCoreData &densestCore) {
  graph->maxCliquecore = 0;
  double max = coreResults[0][2];

  // Parallelize max finding
  double localMax = coreResults[0][2];
  int localMaxCliqueCore = 0;

#pragma omp parallel for reduction(max : localMax, localMaxCliqueCore)
  for (ui i = 1; i < graph->n; i++) {
    if (localMax < coreResults[i][2]) {
      localMax = coreResults[i][2];
    }
    if (localMaxCliqueCore < coreResults[i][1]) {
      localMaxCliqueCore = coreResults[i][1];
    }
  }

  max = localMax;
  graph->maxCliquecore = localMaxCliqueCore;

  int lowerBound = (int)ceil(max);
  int index = 1;
  vector<int> deletedVertices;
  deletedVertices.resize(graph->n, 0);

  // This loop is sequential due to dependency on index
  for (; index < graph->n; index++) {
    if (coreResults[index][1] >= lowerBound) {
      deletedVertices.push_back(coreResults[index][0]);
      break;
    }
    deletedVertices[(int)coreResults[index][0]] = -1;
  }

  int temp = 0;
  // This loop has dependencies, cannot parallelize easily
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
  // Sequential due to temp counter
  for (ui i = index; i < graph->n; i++) {
    int m = (int)coreResults[i][0];
    newCoreResults[temp][0] = deletedVertices[m];
    newCoreResults[temp][1] = coreResults[i][1];
    temp++;
  }

  // Parallelize graph construction
#pragma omp parallel for schedule(dynamic)
  for (ui i = 0; i < graph->n; i++) {
    if (deletedVertices[i] != -1) {
      for (ui j = 0; j < graph->adjacencyList[i].size(); j++) {
        int neighbor = graph->adjacencyList[i][j];
        if (deletedVertices[neighbor] != -1) {
#pragma omp critical
          {
            densestCore.graph[deletedVertices[i]].push_back(
                deletedVertices[neighbor]);
          }
        }
      }
    }
  }

  densestCore.reverseMap.resize(newGraphSize, 0);

  // Parallelize reverse map construction
#pragma omp parallel for
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
    vector<ui> &cliqueDegree, ui motifSize) {
  vector<ui> reverseCoreSortedVertices(newGraph.size());
  vector<ui> order;
  reverseCoreSortedVertices.resize(newGraph.size(), 0);
  vector<ui> degree;
  degree.resize(newGraph.size(), 0);
  for (ui i = 0; i < newGraph.size(); i++) {
    degree[i] = newGraph[i].size();
  }
  vector<ui> core;
  coreDecompose(newGraph, reverseCoreSortedVertices, degree, core, false);

  order.resize(newGraph.size(), 0);
  for (ui i = 0; i < reverseCoreSortedVertices.size(); i++) {
    order[reverseCoreSortedVertices[i]] = i + 1;
  }
  vector<vector<ui>> DAG;
  generateDAG(newGraph, DAG, order);
  vector<ui> partialClique;
  vector<ui> candidates;
  for (ui i = 0; i < newGraph.size(); i++) {
    candidates.push_back(i);
  }
  vector<ui> label;
  label.resize(newGraph.size(), motif->size);
  vector<ui> validNeighborCount;
  validNeighborCount.resize(newGraph.size(), 0);
  listCliqueRecord(motif->size, partialClique, candidates, label, DAG,
                   validNeighborCount, cliqueData, cliqueDegree);
}

void CDS::listCliqueRecord(ui k, vector<ui> &partialClique,
                           vector<ui> &candidates, vector<ui> &label,
                           vector<vector<ui>> &DAG,
                           vector<ui> &validNeighborCount,
                           unordered_map<string, vector<int>> &cliqueData,
                           vector<ui> cliqueDegree) {
  if (k == 2) {
    string cliqueString = "";
    for (ui i = 0; i < partialClique.size(); i++) {
      cliqueString += to_string(partialClique[0]) + " ";
    }

    int cliqueCount = 0;
    for (ui i = 0; i < candidates.size(); i++) {
      int temp = candidates[i];
      for (int j = 0; j < validNeighborCount[temp]; j++) {
        cliqueString = cliqueString + to_string(temp) + to_string(DAG[temp][j]);
        cliqueCount++;
        cliqueDegree[DAG[temp][j]]++;
        cliqueDegree[temp]++;
        vector<int> tempArr(motif->size + 1);
        for (ui x = 0; x < partialClique.size(); x++) {
          tempArr[x] = partialClique[x];
        }
        tempArr[motif->size - 2] = temp;
        tempArr[motif->size - 1] = DAG[temp][j];
        tempArr[motif->size] = 1;

        cliqueData[cliqueString] = tempArr;
      }
    }

    for (ui i = 0; i < candidates.size(); i++) {
      int temp = candidates[i];
      cliqueDegree[temp] += cliqueCount;
    }

  } else {
    for (int i = 0; i < candidates.size(); i++) {
      int temp = candidates[i];
      vector<ui> validNeighbors;
      for (int j = 0; j < DAG[temp].size(); j++) {
        if (label[DAG[temp][j]] == k) {
          label[DAG[temp][j]] = k - 1;
          validNeighbors.push_back(DAG[temp][j]);
        }
      }

      for (int j = 0; j < validNeighbors.size(); j++) {

        ui canTemp = validNeighbors[j];
        int index = 0;
        for (ui m = DAG[canTemp].size() - 1; m > index; --m) {
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

        if (DAG[canTemp].size() != 0 && label[DAG[canTemp][index]] == k - 1)
          index++;

        validNeighborCount[canTemp] = index;
      }
      partialClique.push_back(temp);
      listCliqueRecord(k - 1, partialClique, candidates, label, DAG,
                       validNeighborCount, cliqueData, cliqueDegree);
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

  vector<unordered_map<int, int>> validEdges(oldGraph.size());

  // Build valid edges set from clique data
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

  newGraph.resize(oldGraph.size());
  int totalEdges = 0;

  // Parallelize edge filtering - each vertex processed independently
#pragma omp parallel for reduction(+ : totalEdges) schedule(dynamic)
  for (ui i = 0; i < oldGraph.size(); i++) {
    for (ui j = 0; j < oldGraph[i].size(); j++) {
      if (validEdges[i].find(oldGraph[i][j]) != validEdges[i].end()) {
        newGraph[i].push_back(oldGraph[i][j]);
        totalEdges += 1;
      }
    }
  }

  return totalEdges / 2;
}
void CDS::BFS(vector<ui> status, int vertex, int index,
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

  vector<ui> status;
  status.resize(newGraph.size(), 0);
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
      conComp.totalCliques = temp[temp.size() - 1];
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

    for (ui i = 1; i <= index; i++) {
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
    for (ui i = 1; i < index + 1; i++) {
      graphList[i].resize(graphSizeList[i]);
    }

    for (ui i = 0; i < newGraph.size(); i++) {
      for (int j = 0; j < newGraph[i].size(); j++) {
        graphList[status[i]][oldToNew[i]].push_back(oldToNew[newGraph[i][j]]);
      }
    }

    conCompList.resize(index);
#pragma omp parallel for schedule(dynamic)
    for (ui i = 1; i <= index; i++) {
      ConnectedComponentData conComp;
      conComp.totalCliques = 0;
      conComp.cliqueDegree.resize(graphSizeList[i], 0);

      for (const auto &entry : cliqueDataList[i]) {
        const vector<int> &temp = entry.second;
        for (ui j = 0; j < temp.size() - 1; j++) {
          conComp.cliqueDegree[temp[j]] += temp[temp.size() - 1];
        }
        conComp.totalCliques = temp[temp.size() - 1];
      }

      conComp.graph = graphList[i];
      conComp.size = graphSizeList[i];
      conComp.cliqueData = cliqueDataList[i];
      conComp.density = (double)conComp.totalCliques / (double)conComp.size;
      conComp.reverseMap = reverseMapList[i];

      // SAFE: unique index per thread
      conCompList[i - 1] = std::move(conComp);
    }
  }
}

void CDS::dynamicExact(vector<ConnectedComponentData> &conCompList,
                       DensestCoreData &densestCore,
                       finalResult &densestSubgraph, bool ub1, bool ub2) {

  // ---- Initial lower bound from components ----
  double globalLowerBound = 0.0;
  ConnectedComponentData baseComp = conCompList[0];

  for (ui i = 0; i < conCompList.size(); i++) {
    if (conCompList[i].density > globalLowerBound) {
      globalLowerBound = conCompList[i].density;
      baseComp = conCompList[i];
    }
  }

  if (ceil(globalLowerBound) < ceil(densestCore.density)) {
    globalLowerBound = densestCore.density;
  }

  double globalUpperBound = densestCore.maxCliqueCore;

  // ---- Initialize global best result ----
  finalResult bestGlobal;
  bestGlobal.density = densestSubgraph.density;
  bestGlobal.size = densestSubgraph.size;
  bestGlobal.verticies = densestSubgraph.verticies;

  // ---- Shared variable for dynamic lower bound updates ----
  double sharedLowerBound = globalLowerBound;
  omp_lock_t lowerBoundLock;
  omp_init_lock(&lowerBoundLock);

  // ---- Parallel region ----
#pragma omp parallel
  {
    finalResult localBest;
    localBest.density = -1.0;

#pragma omp for schedule(dynamic)
    for (ui idx = 0; idx < conCompList.size(); idx++) {

      ConnectedComponentData current = conCompList[idx];

      // Read the current shared lower bound
      double lowerBound;
      omp_set_lock(&lowerBoundLock);
      lowerBound = sharedLowerBound;
      omp_unset_lock(&lowerBoundLock);

      double upperBound = globalUpperBound;

      // ---------- UB1 ----------
      if (ub1) {
        double k = (double)motif->size;
        double log_fact = 0.0;
        for (ui j = 1; j <= motif->size; ++j)
          log_fact += log((double)j);

        double dem = exp(log_fact / k);
        double num = pow((double)current.totalCliques, (k - 1.0) / k);

        upperBound = min(upperBound, num / dem);
      }

      // ---------- UB2 ----------
      if (ub2) {
        double dem = (double)(motif->size * current.totalCliques);
        double ub2val = 0.0;

        for (ui deg = 0; deg < current.cliqueDegree.size(); deg++) {
          double pv = (double)current.cliqueDegree[deg] / dem;
          ub2val += pv * pv;
        }
        upperBound = min(upperBound, ub2val);
      }

      // ---------- Early pruning based on upper bound ----------
      // If the upper bound is less than or equal to current lower bound,
      // this component cannot improve the solution
      if (upperBound <= lowerBound) {
        continue; // Skip this component
      }

      // ---------- Exact solver ----------
      vector<int> res;
      exact(res, current, densestCore, densestSubgraph, upperBound, lowerBound);

      // ---------- Compute density ----------
      long cliqueCount = 0;
      ui vertexCount = 0;

      for (const auto &entry : current.cliqueData) {
        const vector<int> &temp = entry.second;
        bool valid = true;

        for (ui j = 0; j < temp.size() - 1; j++) {
          if (res[temp[j]] == -1) {
            valid = false;
            break;
          }
        }

        if (valid)
          cliqueCount += temp.back();
      }

      for (ui j = 0; j < current.size; j++) {
        if (res[j] != -1)
          vertexCount++;
      }

      if (vertexCount == 0)
        vertexCount = current.size;

      double density = (double)cliqueCount / (double)vertexCount;

      // ---------- Update local best ----------
      if (density > localBest.density) {
        localBest.density = density;
        localBest.size = vertexCount;
        localBest.verticies.clear();

        for (ui j = 0; j < res.size(); j++) {
          if (res[j] != -1) {
            localBest.verticies.push_back(
                densestCore.reverseMap[current.reverseMap[j]]);
          }
        }

        // ---------- Update shared lower bound if we found a better solution
        // ----------
        omp_set_lock(&lowerBoundLock);
        if (density > sharedLowerBound) {
          sharedLowerBound = density;
        }
        omp_unset_lock(&lowerBoundLock);
      }
    }

    // ---------- Reduce into global best ----------
#pragma omp critical
    {
      if (localBest.density > bestGlobal.density) {
        bestGlobal = localBest;
      }
    }
  }

  // Clean up the lock
  omp_destroy_lock(&lowerBoundLock);

  // ---- Final result ----
  densestSubgraph = bestGlobal;
}
void CDS::exact(vector<int> res, ConnectedComponentData &conComp,
                DensestCoreData &densestCore, finalResult &densestSubgraph,
                double upperBound, double lowerBound) {

  double alpha = (upperBound + lowerBound) / 2;
  double bais = 1.0 / (conComp.size * (conComp.size - 1));
  if (bais < 0.000000000000001) {
    bais = 0.000000000000001;
  }

  vector<unordered_map<int, array<double, 2>>> flowNetwork;
  vector<int> parent;

  createFlownetwork(flowNetwork, conComp, alpha);

  res.clear();
  res.resize(conComp.size, 1);
  while ((upperBound - lowerBound) > bais) {
    double currentDesnisty = edmondsKarp(flowNetwork, parent, conComp, alpha);

    if (currentDesnisty == conComp.totalCliques * motif->size) {
      upperBound = alpha;
    } else {
      lowerBound = alpha;
      for (ui i = 0; i < conComp.size; i++) {
        res[i] = parent[i];
      }
    }
    alpha = (upperBound + lowerBound) / 2;
    updateFlownetwork(flowNetwork, conComp, alpha);
  }
}

void CDS::createFlownetwork(
    vector<unordered_map<int, array<double, 2>>> &flowNetwork,
    ConnectedComponentData &conComp, double alpha) {
  int flowNetworkSize = conComp.size + conComp.totalCliques + 2;
  int a = conComp.totalCliques;
  flowNetwork.clear();
  flowNetwork.resize(flowNetworkSize);
  int i = 0;
  double weight = 0.0;
  i = conComp.size;
  for (const auto &entry : conComp.cliqueData) {
    const vector<int> &temp = entry.second;
    weight = temp[motif->size] * (motif->size - 1);
    for (ui x = 0; x < motif->size; x++) {
      array<double, 2> temp1 = {static_cast<double>(temp[motif->size]),
                                static_cast<double>(temp[motif->size])};
      array<double, 2> temp2 = {weight, weight};

      // add edge from motif to vertex
      flowNetwork[i][temp[a]] = temp2;

      // add edge from vertex to motif
      flowNetwork[temp[a]][i] = temp1;
    }
    ++i;
  }

  int source = conComp.totalCliques + conComp.size;
  int sink = source + 1;
  for (i = 0; i < conComp.size; i++) {
    array<double, 2> temp1 = {0.0, 0.0};
    flowNetwork[i][source] = temp1;
    array<double, 2> temp2 = {alpha * (motif->size), alpha * (motif->size)};
    flowNetwork[i][sink] = temp2;
    array<double, 2> temp3 = {static_cast<double>(conComp.cliqueDegree[i]),
                              static_cast<double>(conComp.cliqueDegree[i])};
    flowNetwork[source][i] = temp3;
    array<double, 2> temp4 = {0.0, 0.0};
    flowNetwork[sink][i] = temp4;
  }
}

void CDS::updateFlownetwork(
    vector<unordered_map<int, array<double, 2>>> &flowNetwork,
    ConnectedComponentData &conComp, double alpha) {
  int sink = conComp.size + conComp.totalCliques + 1;
  // reset available capacity = max capacity for all edges
  for (int i = 0; i <= sink; ++i) {
    for (auto &entry : flowNetwork[i]) {
      auto &temp_array = entry.second;
      temp_array[0] = temp_array[1];
    }
  }

  // update edges from graph vertices to sink
  for (int i = 0; i < conComp.size; ++i) {
    auto &temp_array = flowNetwork[i][sink];
    temp_array[0] = alpha * motif->size;
    temp_array[1] = alpha * motif->size;
  }
}

double
CDS::edmondsKarp(vector<unordered_map<int, array<double, 2>>> &flowNetwork,
                 vector<int> &parent, ConnectedComponentData &conComp,
                 double alpha) {
  parent.clear();
  parent.resize(flowNetwork.size());
  double minCut = augmentPath(flowNetwork, parent, conComp, alpha);

  double sum = 0;
  vector<double> temp;
  int sink = conComp.size + conComp.totalCliques + 1;
  int source = sink - 1;
  while (minCut != -1) {
    int cur = sink;
    while (cur != source) {
      flowNetwork[parent[cur]][cur][0] =
          flowNetwork[parent[cur]][cur][0] - minCut;
      flowNetwork[cur][parent[cur]][0] =
          flowNetwork[cur][parent[cur]][0] + minCut;
      cur = parent[cur];
    }
    sum += minCut;
    minCut = augmentPath(flowNetwork, parent, conComp, alpha);
  }
  return sum;
}

double
CDS::augmentPath(vector<unordered_map<int, array<double, 2>>> &flowNetwork,
                 vector<int> &parent, ConnectedComponentData &conComp,
                 double alpha) {
  double maxflow = DINF;
  fill(parent.begin(), parent.end(), -1);
  int source = conComp.size + conComp.totalCliques;
  int sink = source + 1;
  queue<int> queue;
  queue.push(source);
  parent[source] = source;
  while (!queue.empty()) {
    int p = queue.front();
    queue.pop();
    if (p == sink) {
      while (p != source) {
        if (maxflow > flowNetwork[parent[p]][p][0]) {
          maxflow = flowNetwork[parent[p]][p][0];
        }
        p = parent[p];
      }
      break;
    }
    for (auto entry : flowNetwork[p]) {
      if (parent[entry.first] == -1 && entry.second[0] > 0) {
        parent[entry.first] = p;
        queue.push(entry.first);
      }
    }
  }

  if (parent[sink] == -1) {
    return -1;
  }
  return maxflow;
}

void CDS::DSD() {
  vector<vector<double>> results;
  cliqueCoreDecompose(results);
  if (debug) {
    for (ui i = 0; i < results.size(); i++) {
      cout << "Vertex: " << results[i][0] << " Clique Degree: " << results[i][1]
           << " Density: " << results[i][2]
           << " Total Cliques remaining: " << results[i][3]
           << " Clique Core value: " << results[i][4] << endl;
    }
  }

  DensestCoreData densestCore;
  locateDensestCore(results, densestCore);
  unordered_map<string, vector<int>> cliqueData;
  vector<ui> cliqueDegree;
  cliqueEnumerationListRecord(densestCore.graph, cliqueData, cliqueDegree,
                              motif->size);

  vector<vector<ui>> newGraph;

  int validEdgeCount =
      pruneInvalidEdges(densestCore.graph, newGraph, cliqueData);

  vector<ConnectedComponentData> conCompList;
  connectedComponentDecompose(newGraph, cliqueData, conCompList);
  finalResult densestSubgraph;
  dynamicExact(conCompList, densestCore, densestSubgraph, false, false);
  cout << "Results" << endl;
  cout << "K: " << motif->size << endl;
  cout << "Density: " << densestSubgraph.density << endl;
  cout << "Size: " << densestSubgraph.size << endl;
}