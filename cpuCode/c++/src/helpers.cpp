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
  vector<ui> tempArray(graph->n);
  iota(order.begin(), order.end(), 0);

  // sort vertices by core value descending
  sort(tempArray.begin(), tempArray.end(),
       [&](int a, int b) { return graph->core[a] > graph->core[b]; });

  order.resize(graph->n, 0);
  for (ui i = 0; i < tempArray.size(); i++) {
    order[i] = i + 1;
  }
}