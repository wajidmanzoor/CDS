#pragma once

#include "graph.h"

class CDS {
private:
  Graph *graph;
  Motif *motif;
  bool ub1;
  bool ub2;
  void get2Dneighborhood(unordered_map<int, long> &subgraphResults, int index,
                         vector<int> &mark, vector<int> &array,
                         vector<int> &map_s);
  void cliqueEnumerationFast();
  void getlistingOrder(vector<ui> &order);
  void generateDAG(const vector<vector<ui>> adjList, vector<vector<ui>> &DAG,
                   vector<ui> &order);
  void coreDecompose(const vector<vector<ui>> adjList,
                     vector<ui> &reverseCoreSortedVertices, vector<ui> &degree,
                     vector<ui> &core, bool fullGraph = true);
  void listCliques(ui k, vector<ui> &partialClique, vector<ui> &candidates,
                   vector<ui> &label, vector<vector<ui>> &DAG,
                   vector<ui> &validNeighborCount);
  void cliqueEnumerationSubgraph(vector<vector<ui>> &subGraph,
                                 vector<ui> &subGraphCliqueDegree, ui motifSize,
                                 ui vertex);
  void listCliqueContainsVertex(ui k, vector<ui> &partialClique,
                                vector<ui> &candidates, vector<ui> &label,
                                vector<vector<ui>> &DAG,
                                vector<ui> &validNeighborCount,
                                vector<ui> &cliqueDegree, ui vertex);
  void listCliqueRecord(ui k, vector<ui> &partialClique, vector<ui> &candidates,
                        vector<ui> &label, vector<vector<ui>> &DAG,
                        vector<ui> &validNeighborCount,
                        unordered_map<string, vector<int>> &cliqueData,
                        vector<long> cliqueDegree);
  void BFS(vector<ui> &status, int vertex, int index,
           const vector<vector<ui>> &newGraph);
  void exact(vector<int> &res, ConnectedComponentData &C, float lowerBound,
             float upperBound);
  void createFlownetwork(FlowNetwork &FN, ConnectedComponentData &C,
                         float alpha);
  // float edmondsKarp(FlowNetwork &FN, vector<int> &parent);
  float edmondsKarp(FlowNetwork &FN, vector<pair<int, int>> &parent);

  void updateFlownetwork(FlowNetwork &FN, ConnectedComponentData &C,
                         float alpha);

  // float augmentPath(FlowNetwork &FN, vector<int> &parent);
  float augmentPath(FlowNetwork &FN, vector<pair<int, int>> &parent);

public:
  CDS();
  CDS(Graph *graph, Motif *motif, bool ub1, bool ub2);
  void cliqueCoreDecompose(vector<vector<double>> &results);
  void
  cliqueEnumerationListRecord(vector<vector<ui>> newGraph,
                              unordered_map<string, vector<int>> &cliqueData,
                              vector<long> &cliqueDegree, ui motifSize);
  void locateDensestCore(vector<vector<double>> &coreResults,
                         DensestCoreData &densestCore);

  int pruneInvalidEdges(vector<vector<ui>> &oldGraph,
                        vector<vector<ui>> &newGraph,
                        unordered_map<string, vector<int>> &cliqueData);

  void
  connectedComponentDecompose(vector<vector<ui>> &newGraph,
                              unordered_map<string, vector<int>> &cliqueData,
                              vector<ConnectedComponentData> &conCompList);
  void dynamicExact(vector<ConnectedComponentData> &conCompList,
                    DensestCoreData &densestCore, finalResult &densestSubgraph,
                    bool ub1, bool ub2);

  void DSD();
};