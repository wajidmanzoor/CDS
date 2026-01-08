#pragma once

#include "graph.h"

class CDS {
private:
  Graph *graph;
  Motif *motif;
  unordered_map<int, long> get2Dneighborhood(int index, vector<int> mark,
                                             vector<int> array,
                                             vector<int> map_s);
  void cliqueEnumerationFast();
  void getlistingOrder(vector<ui> &order);
  void generateDAG(vector<vector<ui>> &DAG, vector<ui> &order);
  void coreDecompose(vector<ui> &reverseCoreSortedVertices);
  void listCliques(ui k, vector<ui> &partialClique, vector<ui> &candidates,
                   vector<ui> &label, vector<vector<ui>> DAG,
                   vector<ui> validNeighborCount);

public:
  CDS();
  CDS(Graph *graph, Motif *motif);
  vector<vector<double>> cliqueCoreDecompose();
  void cliqueEnumerationListRecord();
};