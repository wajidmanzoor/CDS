#pragma once

#include "graph.h"

class CDS {
private:
  Graph *graph;
  Motif *motif;
  void get2Dneighborhood(unordered_map<int, long> &subgraphResults, int index,
                         vector<int> mark, vector<int> array,
                         vector<int> map_s);
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
  void CDS::listCliqueContainsVertex(ui k, vector<ui> &partialClique,
                                     vector<ui> &candidates, vector<ui> &label,
                                     vector<vector<ui>> &DAG,
                                     vector<ui> &validNeighborCount,
                                     vector<ui> &cliqueDegree, ui vertex);
  void listCliqueRecord(ui k, vector<ui> &partialClique, vector<ui> &candidates,
                        vector<ui> &label, vector<vector<ui>> &DAG,
                        vector<ui> &validNeighborCount,
                        unordered_map<string, vector<int>> &cliqueData,
                        vector<ui> cliqueDegree);

public:
  CDS();
  CDS(Graph *graph, Motif *motif);
  void cliqueCoreDecompose(vector<vector<double>> &results);
  void
  cliqueEnumerationListRecord(vector<vector<ui>> newGraph,
                              unordered_map<string, vector<int>> &cliqueData,
                              vector<ui> &cliqueDegree, ui motifSize);
  void locateDensestCore(vector<vector<double>> &coreResults,
                         DensestCoreData &densestCore);

  int pruneInvalidEdges(vector<vector<ui>> &oldGraph,
                        vector<vector<ui>> &newGraph,
                        unordered_map<string, vector<int>> &cliqueData);

  void DSD();
};