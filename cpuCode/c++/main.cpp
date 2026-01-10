#pragma once
#include "inc/common.h"
#include "inc/graph.h"
#include "inc/helpers.h"
int main() {

  Graph graph("data/test_graph.txt");
  Motif motif("/data/test_motif.txt");

  CDS cds(&graph, &motif);
  cds.DSD();
}