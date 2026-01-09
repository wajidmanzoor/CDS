#pragma once
#include "../inc/helpers.h"
#include "../inc/graph.h"
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
  for (ui i = 0; i < graph->n; i++) {
    graph->cliqueCore[i] = graph->cliqueDegree[i];
  }

  int totalCliques = 0;

  graph->maxCliqueDegree = 0;
  for (ui i = 0; i < graph->n; i++) {
    if (graph->cliqueDegree[i] > graph->maxCliqueDegree) {
      graph->maxCliqueDegree = graph->cliqueDegree[i];
    }
    totalCliques += graph->cliqueDegree[i];
  }

  totalCliques = totalCliques / static_cast<double>(motif->size);

  // data structure used to save clique core decompose results

  results.resize(graph->n + 1, vector<double>(5, 0.0));
  results[0][2] = totalCliques / static_cast<double>(graph->n);
  results[0][3] = totalCliques;
  // Create bins for counting the number of vertices with each clique degree
  vector<long> bins;
  bins.resize(graph->maxCliqueDegree + 1, 0);

  for (ui i = 0; i < graph->n; i++) {
    bins[graph->cliqueDegree[i]]++;
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
    int index = 0;
    long indexMin = 0xFFFFFF;
    for (ui j = 0; j < graph->n; j++) {
      if (indexMin > graph->cliqueDegree[j] && mark[j] == 0) {
        indexMin = graph->cliqueDegree[j];
        index = j;
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
                            int index, vector<int> mark, vector<int> array,
                            vector<int> map_s) {
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
  for (ui i = 0; i < adjList.size(); i++) {
    count = 0;
    for (ui j = 0; j < adjList[i].size(); j++) {
      if (order[adjList[i][j]] > order[i]) {
        count++;
      }
    }
    DAG[i].resize(count, 0);
    int index = 0;
    for (ui j = 0; j < adjList[i].size(); j++) {
      if (order[adjList[i][j]] > order[i]) {
        DAG[i][index] = adjList[i][j];
        index++;
      }
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
  for (ui i = 1; i < graph->n; i++) {
    if (max < coreResults[i][2]) {
      max = coreResults[i][2];
    }
    if (graph->maxCliquecore < coreResults[i][1]) {
      graph->maxCliquecore = coreResults[i][1];
    }
  }

  int lowerBound = (int)ceil(max);

  int index = 1;
  vector<int> deletedVertices;
  deletedVertices.resize(graph->n, 0);
  for (; index < graph->n; index++) {
    if (coreResults[index][1] >= lowerBound) {
      deletedVertices.push_back(coreResults[index][0]);
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

int pruneInvalidEdges(vector<vector<ui>> &oldGraph,
                      vector<vector<ui>> &newGraph,
                      unordered_map<string, vector<int>> &cliqueData) {
  int count = 0;
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
    vector<vector<ui>> &newGraph, vector<ConnectedComponentData> &conComps) {

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

  } else {
  }
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
}