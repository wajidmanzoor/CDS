#pragma once

#include "graph.h"

class CDS {
private:
  Graph *graph;
  Motif *motif;
  unordered_map<int, long> get2Dneighborhood(int index, vector<int> mark,
                                             vector<int> array,
                                             vector<int> map_s);
  vector<long> cliqueEnumerationFast();
  void getlistingOrder(vector<ui> &order);
  vector<vector<ui>> generateDAG();

public:
  CDS();
  CDS(Graph *graph, Motif *motif);
  vector<vector<double>> cliqueCoreDecompose();
  void cliqueEnumerationListRecord();
};