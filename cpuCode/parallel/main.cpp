#include "inc/common.h"
#include "inc/graph.h"
#include "inc/helpers.h"
int main() {

  omp_set_num_threads(THREAD_COUNT);

  Graph graph("data/test_graph.txt");
  Motif motif("/data/test_motif.txt");

  CDS cds(&graph, &motif);
  cds.DSD();
}