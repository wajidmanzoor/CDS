#include "inc/common.h"
int main() {

  Graph graph("data/test_graph.txt");
  Motif motif("data/test_motif.txt");

  DynamicExactAlgo dynamic(&graph, &motif);
  DynamicExactAlgo.listAllCliques();
  DynamicExactAlgo.cliqueCoreDecompose();
  DynamicExactAlgo.findDensestCore();
}