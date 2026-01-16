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
  newArray.resize(graph->n, 0);
  vector<int> newMap;
  newMap.resize(graph->n, 0);
  // cout << "before clique enum" << endl;

  cliqueEnumerationFast();

  // cout << "after clique enum" << endl;
  unordered_map<int, long> twoDNeighborhood;
  graph->cliqueCore.resize(graph->n, 0);
  for (ui i = 0; i < graph->n; i++) {
    graph->cliqueCore[i] = graph->cliqueDegree[i];
  }
  // cout << "here" << endl;

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

  // cout << "total clique " << totalCliques << endl;

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
    results[count][1] = graph->cliqueDegree[index];
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

        deleteCount = deleteCount / (motif->size);
        if (debug)
          cout << "deleted cliques: " << deleteCount << endl;
        totalCliques -= deleteCount;
        results[count][3] = totalCliques;
        if (debug)
          cout << "deleted verticies: " << count << endl;
        if (graph->n - count > 0) {

          results[count][2] =
              totalCliques / static_cast<double>(graph->n - count);
        } else {
          results[count][2] = 0.0;
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
  for (ui i = 0; i < count; i++) {
    int vertex = tempList[i];
    mapArray[i] = vertex;
    map_s[vertex] = num;
    int tempCount = 0;
    for (int j = 0; j < graph->adjacencyList[vertex].size(); j++) {
      int neighbor = graph->adjacencyList[vertex][j];
      if (array[neighbor] != 0 && neighbor != vertex) {
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
    // cout << "neigh of " << i << " : ";
    for (int j = 0; j < subGraph[i].size(); ++j) {
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
  graph->cliqueDegree.resize(graph->n, 0);
  graph->totalCliques = 0;

  // cout << "Before listCliques" << endl;

  listCliques(motif->size, partialClique, candidates, label, DAG,
              validNeighborCount);

  // cout << "donefinal" << endl;
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
        // cout << "i: " << i << " j: " << index << " DAG: " << DAG[i][index];
        DAG[i][index] = adjList[i][j];
        // cout << " new DAG " << DAG[i][index] << endl;
        index++;
      }
    }
  }
}

void CDS::listCliques(ui k, vector<ui> &partialClique, vector<ui> &candidates,
                      vector<ui> &label, vector<vector<ui>> &DAG,
                      vector<ui> &validNeighborCount) {

  if (debug) {
    cout << "partial Clique: ";
    for (ui pc : partialClique) {
      cout << pc << " ";
    }
    cout << endl << " candidates: ";
    for (ui c : candidates) {
      cout << c << " ";
    }
    cout << endl << " label: ";
    for (ui l : label) {
      cout << l << " ";
    }
    cout << endl << " vnc: ";
    for (ui vn : validNeighborCount) {
      cout << vn << " ";
    }
    cout << endl;
    int x = 0;
    for (ui vn : validNeighborCount) {
      cout << "valid neighbors of " << x << ": ";
      if (vn > 0) {
        for (ui neig : DAG[x]) {
          cout << neig << " ";
        }
      }
      x++;

      cout << endl;
    }
    cout << "total Cliques " << graph->totalCliques << endl;
  }

  if (k == 2) {
    if (debug)
      cout << "----------------------" << endl;
    string cliqueString = "";
    for (ui i = 0; i < partialClique.size(); i++) {
      cliqueString += to_string(partialClique[i]) + " ";
    }
    // cout << "I am here" << endl;

    long cliqueCount = 0;
    // cout << "can size " << candidates.size() << endl;
    for (ui i = 0; i < candidates.size(); i++) {

      int temp = candidates[i];
      // cout << "k-1 can:" << temp << " vn size" << validNeighborCount[temp]
      //    << endl;
      for (int j = 0; j < validNeighborCount[temp]; j++) {
        string wajid = cliqueString + to_string(candidates[i]) + " " +
                       to_string(DAG[temp][j]);
        if (debug)
          cout << "clique number " << graph->totalCliques << " : " << wajid
               << endl;

        // cout << "j " << endl;
        cliqueCount++;
        graph->totalCliques++;
        graph->cliqueDegree[DAG[temp][j]]++;
        graph->cliqueDegree[temp]++;
      }
    }

    // /cout << cliqueString << " Total Cliques " << graph->totalCliques <<
    // endl;
    for (ui i = 0; i < partialClique.size(); i++) {
      int temp = partialClique[i];
      // cout << " temp " << temp << " clique Count " << cliqueCount << " Prev "
      //     << graph->cliqueDegree[temp] << endl;
      graph->cliqueDegree[temp] += cliqueCount;
    }
    if (debug)
      cout << "----------------------" << endl;

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
      // cout << "Can :" << temp << endl;

      for (int j = 0; j < validNeighbors.size(); j++) {

        ui canTemp = validNeighbors[j];
        // cout << "Valid Neighs: " << canTemp << endl;

        int index = 0;
        for (int m = DAG[canTemp].size() - 1; m > index; --m) {
          if (label[DAG[canTemp][m]] == k - 1) {
            while (index < m && label[DAG[canTemp][index]] == k - 1) {
              // cout << "index: " << index << endl;
              index++;
            }
            if (label[DAG[canTemp][index]] != k - 1) {
              // cout << "final index :" << index << " m " << m << endl;
              int temp1 = DAG[canTemp][m];
              DAG[canTemp][m] = DAG[canTemp][index];
              DAG[canTemp][index] = temp1;
            }
          }
        }

        // cout << "index at end: " << index << endl;

        // cout << "DAG of " << canTemp << " size: " << DAG[canTemp].size()
        // << endl;

        if (DAG[canTemp].size() != 0)
          if (label[DAG[canTemp][index]] == k - 1)
            index++;
        // cout << "index after: " << index << " label " << canTemp << endl;

        validNeighborCount[canTemp] = index;
      }

      // cout << "done" << endl;
      partialClique.push_back(temp);
      listCliques(k - 1, partialClique, validNeighbors, label, DAG,
                  validNeighborCount);
      partialClique.pop_back();
      // cout << "done2" << endl;
      // cout << " size vn: " << validNeighbors.size() << endl;
      for (ui j = 0; j < validNeighbors.size(); j++) {
        // cout << "j " << j << endl;
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
  // cout << "before clique list subgraph clique degree " << endl;
  /*for (ui v : subGraphCliqueDegree) {
    cout << v << " ";
  }
  cout << endl;*/

  listCliqueContainsVertex(motif->size, partialClique, candidates, label, DAG,
                           validNeighborCount, subGraphCliqueDegree, 0);
}

void CDS::listCliqueContainsVertex(ui k, vector<ui> &partialClique,
                                   vector<ui> &candidates, vector<ui> &label,
                                   vector<vector<ui>> &DAG,
                                   vector<ui> &validNeighborCount,

                                   vector<ui> &cliqueDegree, ui vertex) {
  /*if (debug) {
    cout << "partial Clique: ";
    for (ui pc : partialClique) {
      cout << pc << " ";
    }
    cout << endl << " candidates: ";
    for (ui c : candidates) {
      cout << c << " ";
    }
    cout << endl << " label: ";
    for (ui l : label) {
      cout << l << " ";
    }
    cout << endl << " vnc: ";
    for (ui vn : validNeighborCount) {
      cout << vn << " ";
    }
    cout << endl;
    int x = 0;
    for (ui vn : validNeighborCount) {
      cout << "valid neighbors of " << x << ": ";
      if (vn > 0) {
        for (ui neig : DAG[x]) {
          cout << neig << " ";
        }
      }
      x++;

      cout << endl;
    }
    // cout << "total Cliques " << graph->totalCliques << endl;
  }*/

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

      int temp = candidates[i];
      for (int j = 0; j < validNeighborCount[temp]; j++) {
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
    for (int i = 0; i < candidates.size(); i++) {
      int temp = candidates[i];
      vector<ui> validNeighbors;
      // cout << " intial start " << i << " can " << temp << endl;
      for (int j = 0; j < DAG[temp].size(); j++) {
        if (label[DAG[temp][j]] == k) {
          label[DAG[temp][j]] = k - 1;
          validNeighbors.push_back(DAG[temp][j]);
        }
      }
      // cout << " valid neighs size of " << temp << " is "
      //    << validNeighbors.size() << endl;
      for (int j = 0; j < validNeighbors.size(); j++) {

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

  // cout << "max: " << max << endl;

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
    vector<long> &cliqueDegree, ui motifSize) {
  vector<ui> reverseCoreSortedVertices(newGraph.size());
  reverseCoreSortedVertices.resize(newGraph.size(), 0);

  vector<ui> degree;
  degree.resize(newGraph.size(), 0);
  for (ui i = 0; i < newGraph.size(); i++) {
    degree[i] = newGraph[i].size();
  }
  vector<ui> core;
  coreDecompose(newGraph, reverseCoreSortedVertices, degree, core, false);

  vector<ui> order;
  order.resize(newGraph.size(), 0);
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

  listCliqueRecord(motif->size, partialClique, candidates, label, DAG,
                   validNeighborCount, cliqueData, cliqueDegree);
}

void CDS::listCliqueRecord(ui k, vector<ui> &partialClique,
                           vector<ui> &candidates, vector<ui> &label,
                           vector<vector<ui>> &DAG,
                           vector<ui> &validNeighborCount,
                           unordered_map<string, vector<int>> &cliqueData,
                           vector<long> cliqueDegree) {
  if (debug) {
    cout << "partial Clique: ";
    for (ui pc : partialClique) {
      cout << pc << " ";
    }
    cout << endl << " candidates: ";
    for (ui c : candidates) {
      cout << c << " ";
    }
    cout << endl << " label: ";
    for (ui l : label) {
      cout << l << " ";
    }
    cout << endl << " vnc: ";
    for (ui vn : validNeighborCount) {
      cout << vn << " ";
    }
    cout << endl;
    int x = 0;
    for (ui vn : validNeighborCount) {
      cout << "valid neighbors of " << x << ": ";
      if (vn > 0) {
        for (ui neig : DAG[x]) {
          cout << neig << " ";
        }
      }
      x++;

      cout << endl;
    }
    cout << "total Cliques " << graph->totalCliques << endl;
  }

  if (k == 2) {
    if (debug)
      cout << "----------------------" << endl;
    string cliqueString = "";
    for (ui i = 0; i < partialClique.size(); i++) {
      cliqueString += to_string(partialClique[i]) + " ";
    }

    long cliqueCount = 0;
    for (ui i = 0; i < candidates.size(); i++) {
      int temp = candidates[i];
      for (int j = 0; j < validNeighborCount[temp]; j++) {
        string cliqueString1 =
            cliqueString + to_string(temp) + " " + to_string(DAG[temp][j]);
        if (debug)
          cout << "clique number " << graph->totalCliques << " : "
               << cliqueString1 << endl;
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

        cliqueData[cliqueString1] = tempArr;
      }
    }

    for (ui i = 0; i < partialClique.size(); i++) {
      int temp = partialClique[i];
      cliqueDegree[temp] += cliqueCount;
    }
    if (debug)
      cout << "----------------------" << endl;

  } else {
    for (int i = 0; i < candidates.size(); i++) {
      int temp = candidates[i];
      // cout << "temp " << temp << endl;
      vector<ui> validNeighbors;
      for (int j = 0; j < DAG[temp].size(); j++) {
        if (label[DAG[temp][j]] == k) {
          label[DAG[temp][j]] = k - 1;
          validNeighbors.push_back(DAG[temp][j]);
        }
      }
      // cout << validNeighbors.size() << " vn size" << endl;

      for (int j = 0; j < validNeighbors.size(); j++) {

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
      listCliqueRecord(k - 1, partialClique, validNeighbors, label, DAG,
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

  vector<ui> status;
  status.resize(newGraph.size(), 0);
  int index = 0;
  for (ui i = 0; i < newGraph.size(); i++) {
    if (status[i] == 0) {
      index++;
      BFS(status, i, index, newGraph);
    }
  }
  // cout << "here " << endl;

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

    for (ui i = 1; i < index + 1; i++) {
      ConnectedComponentData conComp;
      conComp.totalCliques = 0;
      conComp.cliqueDegree.resize(graphSizeList[i], 0);
      for (const auto &entry : cliqueDataList[i]) {
        const vector<int> &temp = entry.second;
        for (ui i = 0; i < temp.size() - 1; i++) {
          conComp.cliqueDegree[temp[i]] += temp[temp.size() - 1];
        }
        conComp.totalCliques = temp[temp.size() - 1];
      }
      conComp.graph = graphList[i];
      conComp.size = graphSizeList[i];
      conComp.cliqueData = cliqueDataList[i];
      conComp.density = (double)conComp.totalCliques / ((double)conComp.size);
      conComp.reverseMap = reverseMapList[i];
      conCompList.push_back(conComp);
    }
  }
}

void CDS::dynamicExact(vector<ConnectedComponentData> &conCompList,
                       DensestCoreData &densestCore,
                       finalResult &densestSubgraph, bool ub1, bool ub2) {

  ConnectedComponentData current, index;
  double lowerBound = 0;
  current = conCompList[0];
  for (ui i = 0; i < conCompList.size(); i++) {

    if (lowerBound < conCompList[i].density) {
      lowerBound = current.density;
      current = conCompList[i];
    }
  }
  if (ceil(lowerBound) < ceil(densestCore.density)) {
    lowerBound = densestCore.density;
  }

  double upperBound = densestCore.maxCliqueCore;

  densestSubgraph.verticies.resize(current.size);

  for (ui i = 0; i < current.graph.size(); i++) {
    densestSubgraph.verticies[i] =
        densestCore.reverseMap[current.reverseMap[i]];
  }

  // TODO: add new density bounds

  for (ui i = 0; i < conCompList.size(); i++) {
    current = conCompList[i];

    if (ub1) {
      double k = (double)motif->size;
      double log_fact = 0.0;
      for (ui i = 1; i <= motif->size; ++i)
        log_fact += log((double)i);

      double dem = exp(log_fact / k);
      double num = pow((double)current.totalCliques, (k - 1.0) / k);

      double upperBound1 = num / dem;
      if (upperBound1 < upperBound) {
        upperBound = upperBound1;
      }
    }
    if (ub2) {
      double dem = (double)(motif->size * current.totalCliques);
      double upperBound2 = 0.0;
      for (ui deg = 0; deg < current.cliqueDegree.size(); deg++) {
        double pv = (double)current.cliqueDegree[deg] / dem;
        pv = pv * pv;
        upperBound2 += pv;
      }
      upperBound2 = upperBound2 * current.totalCliques;
      if (upperBound2 < upperBound) {
        upperBound = upperBound2;
      }
    }
    vector<int> res;
    exact(res, current, densestCore, densestSubgraph, upperBound, lowerBound);

    cout << endl;
    long cliqueCount = 0;
    ui vertexCount = 0;
    for (const auto &entry : current.cliqueData) {
      const vector<int> &temp = entry.second;
      // cout << "current clique: ";

      int i = 0;
      for (; i < temp.size() - 1; ++i) {

        if (res[temp[i]] == -1) {
          break;
        }
      }

      if (i == temp.size() - 1) {
        cliqueCount += temp[i];
      }
    }
    // cout << "FINAL CLIQUE COUNT: " << cliqueCount << endl;

    for (ui i = 0; i < current.size; i++) {
      if (res[i] != -1) {
        vertexCount++;
      }
    }
    // cout << "FINAL vertex COUNT: " << vertexCount << endl;

    if (vertexCount == 0) {
      vertexCount = current.size;
    }

    double temp = (double)(cliqueCount / ((double)vertexCount));

    if (temp > densestSubgraph.density) {
      densestSubgraph.density = temp;
      lowerBound = temp;
      densestSubgraph.size = vertexCount;
      densestSubgraph.verticies.clear();
      densestSubgraph.verticies.resize(vertexCount);
      // TODO: CHECK MAYBE WRong
      for (ui i = 0; i < res.size(); i++) {
        if (res[i] != -1) {
          densestSubgraph.verticies[i] =
              densestCore.reverseMap[current.reverseMap[i]];
        }
      }
    }
  }
}

void CDS::exact(vector<int> &res, ConnectedComponentData &conComp,
                DensestCoreData &densestCore, finalResult &densestSubgraph,
                double upperBound, double lowerBound) {

  double alpha = (upperBound + lowerBound) / 2;
  double bais = 1.0 / (conComp.size * (conComp.size - 1));
  if (bais < 0.000000000000001) {
    bais = 0.000000000000001;
  }

  vector<unordered_map<int, array<double, 2>>> flowNetwork;
  vector<int> parent;

  // cout << "here" << endl;
  createFlownetwork(flowNetwork, conComp, alpha);
  // cout << "here 2 " << endl;

  res.clear();
  res.resize(conComp.size, 1);

  while ((upperBound - lowerBound) > bais) {
    double maxflow = edmondsKarp(flowNetwork, parent, conComp, alpha);

    if (abs((conComp.totalCliques * motif->size) - maxflow) < 1e-2) {
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
  // cout << "FN size: " << flowNetworkSize << endl;
  int i = 0;
  double weight = 0.0;
  i = conComp.size;
  for (const auto &entry : conComp.cliqueData) {
    const vector<int> &temp = entry.second;
    weight = temp[motif->size] * (motif->size - 1);
    // cout << "Key: " << entry.first << " Weight " << weight << endl;
    for (ui x = 0; x < motif->size; x++) {
      array<double, 2> temp1 = {static_cast<double>(temp[motif->size]),
                                static_cast<double>(temp[motif->size])};
      array<double, 2> temp2 = {weight, weight};

      // add edge from motif to vertex
      flowNetwork[i][temp[x]] = temp2;

      // add edge from vertex to motif
      flowNetwork[temp[x]][i] = temp1;
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
  // cout << "FN after loop2" << endl;
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
  cliqueEnumerationListRecord(densestCore.graph, cliqueData, cliqueDegree,
                              motif->size);
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

  int validEdgeCount =
      pruneInvalidEdges(densestCore.graph, newGraph, cliqueData);
  if (debug) {
    cout << "Graph after pruning:  " << endl;
    for (ui i = 0; i < newGraph.size(); i++) {
      cout << "Neigh of i: " << i << " : ";
      for (ui j = 0; j < newGraph[i].size(); j++) {
        cout << newGraph[i][j] << " ";
      }
      cout << endl;
    }
  }

  vector<ConnectedComponentData> conCompList;
  connectedComponentDecompose(newGraph, cliqueData, conCompList);
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
  dynamicExact(conCompList, densestCore, densestSubgraph, true, true);
  cout << "Final Results" << endl;
  cout << "K: " << motif->size << endl;
  cout << "Density: " << densestSubgraph.density << endl;
  cout << "Size: " << densestSubgraph.size << endl;
}