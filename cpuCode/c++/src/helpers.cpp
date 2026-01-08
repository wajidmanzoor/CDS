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

vector<vector<double>> CDS ::cliqueCoreDecompose() {

  vector<int> mark;
  mark.resize(graph->n, 0);
  vector<int> arrayIndex;
  arrayIndex.resize(graph->n, 0);
  vector<int> newArray;
  arrayIndex.resize(graph->n, 0);
  vector<int> newMap;
  newMap.resize(graph->n, 0);

  unordered_map<int, long> twoDNeighborhood;
  // TODO add clique enemuration.

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

  vector<vector<double>> results;
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
      twoDNeighborhood = get2Dneighborhood(index, mark, arrayIndex, newMap);
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

  return results;
}

unordered_map<int, long> CDS::get2Dneighborhood(int index, vector<int> mark,
                                                vector<int> array,
                                                vector<int> map_s) {

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

  vector<long> subGraphCliqueDegree;
  subGraphCliqueDegree.resize(count, 0);

  // TODO: add clique enumeration for subgraph and return the count of cliques
  // for each vertex in the subgraph.

  unordered_map<int, long> subgraphResults;
  for (ui i = 0; i < count; i++) {
    if (subGraphCliqueDegree[i] > 0) {
      subgraphResults[mapArray[i]] = subGraphCliqueDegree[i];
    }
  }
  return subgraphResults;
}

void CDS::getlistingOrder(vector<ui> &order) {

  // Todo: add core decompose

  // verticies sorted by reverse core values
  vector<ui> reverseCoreSortedVertices(graph->n);

  coreDecompose(reverseCoreSortedVertices);
  order.resize(graph->n, 0);
  for (ui i = 0; i < reverseCoreSortedVertices.size(); i++) {
    order[reverseCoreSortedVertices[i]] = i + 1;
  }
}

void CDS::coreDecompose(vector<ui> &reverseCoreSortedVertices) {
  reverseCoreSortedVertices.resize(graph->n, 0);
  graph->maxDegree = 0;
  for (ui i = 0; i < graph->n; i++) {
    if (graph->degree[i] > graph->maxDegree) {
      graph->maxDegree = graph->degree[i];
    }
    graph->core[i] = graph->degree[i];
  }

  vector<ui> bins;
  bins.resize(graph->maxDegree + 1, 0);

  for (ui i = 0; i < graph->n; i++) {
    bins[graph->degree[i]]++;
  }

  ui start = 0;
  for (ui i = 0; i <= graph->maxDegree; i++) {
    ui temp = bins[i];
    bins[i] = start;
    start += temp;
  }
  vector<ui> pos;
  pos.resize(graph->n, 0);
  vector<ui> sortedVertices;
  sortedVertices.resize(graph->n, 0);
  for (ui i = 0; i < graph->n; i++) {
    pos[i] = bins[graph->degree[i]];
    sortedVertices[pos[i]] = i;
    bins[graph->degree[i]]++;
  }

  for (ui i = graph->maxDegree; i > 0; i--) {
    bins[i] = bins[i - 1];
  }
  bins[0] = 0;

  for (ui i = 0; i < graph->n; i++) {
    int vertex = sortedVertices[i];
    graph->core[vertex] = graph->degree[vertex];
    for (ui j = 0; j < graph->adjacencyList[vertex].size(); j++) {
      int neighbor = graph->adjacencyList[vertex][j];
      if (graph->core[neighbor] > graph->core[vertex]) {
        int du = graph->core[neighbor];
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
        graph->core[neighbor]--;
      }
      reverseCoreSortedVertices[graph->n - i - 1] = vertex;
    }
  }
}

void CDS::cliqueEnumerationFast() {
  vector<ui> order;
  getlistingOrder(order);
  vector<vector<ui>> DAG;
  generateDAG(DAG, order);
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

void CDS::generateDAG(vector<vector<ui>> &DAG, vector<ui> &order) {
  int count;
  DAG.resize(graph->n);
  for (ui i = 0; i < graph->n; i++) {
    count = 0;
    for (ui j = 0; j < graph->adjacencyList[i].size(); j++) {
      if (order[graph->adjacencyList[i][j]] > order[i]) {
        count++;
      }
    }
    DAG[i].resize(count, 0);
    int index = 0;
    for (ui j = 0; j < graph->adjacencyList[i].size(); j++) {
      if (order[graph->adjacencyList[i][j]] > order[i]) {
        DAG[i][index] = graph->adjacencyList[i][j];
        index++;
      }
    }
  }
}

void CDS::listCliques(ui k, vector<ui> &partialClique, vector<ui> &candidates,
                      vector<ui> &label, vector<vector<ui>> DAG,
                      vector<ui> validNeighborCount) {

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
        graph->cliqueCore[temp]++;
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