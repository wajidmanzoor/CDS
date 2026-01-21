
#include "inc/common.h"
#include "inc/graph.h"
#include "inc/helpers.h"
int main(int argc, const char *argv[]) {
  if (argc != 5) {
    cout << "Server wrong input parameters!" << endl;
    exit(1);
  }
  omp_set_num_threads(THREAD_COUNT);

  string graphpath = argv[1];
  string motifpath = argv[2];
  bool ub1 = atoi(argv[3]) != 0;
  bool ub2 = atoi(argv[4]) != 0;
  Graph graph(graphpath);
  Motif motif(motifpath);

  CDS cds(&graph, &motif, ub1, ub2);
  cds.DSD();
}
