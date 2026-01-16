#include "inc/common.h"
#include "inc/graph.h"
#include "inc/helpers.h"
int main(int argc, const char *argv[]) {
  if (argc != 3) {
    cout << "Server wrong input parameters!" << endl;
    exit(1);
  }
  string graphpath = argv[1];
  string motifpath = argv[2];
  Graph graph(graphpath);
  Motif motif(motifpath);

  CDS cds(&graph, &motif);
  cds.DSD();
}