#include "../inc/graph.h"
#include <numeric>
#include <sstream>

Graph::Graph() {
  // Default constructor implementation
}

Graph::Graph(std::string path) {
  std::string buffer;
  std::ifstream inputFile(path, std::ios::in);

  if (!inputFile.is_open()) {
    std::cout << "Graph file Open Failed " << std::endl;
    exit(1);
  } else {
    std::string line;
    std::getline(inputFile, line);
    std::istringstream iss(line);
    iss >> n >> m;
    adjacencyList.resize(n);

    degree.resize(n, 0);
    int vertex, neigh;
    while (std::getline(inputFile, line)) {
      std::istringstream iss(line);
      iss >> vertex;
      while (iss >> neigh) {
        if (vertex == neigh)
          continue;
        adjacencyList[vertex].push_back(neigh);
      }
      degree[vertex] = adjacencyList[vertex].size();
    }
  }

  inputFile.close();
  if (debug) {
    std::cout << "n=" << n << ", m=" << m << std::endl;

    cout << "adjacencyList " << endl;
    for (ui i = 0; i < n; i++) {
      cout << i << ": ";
      for (ui j = 0; j < adjacencyList[i].size(); j++) {
        cout << adjacencyList[i][j] << " ";
      }
      cout << endl;
    }
  }

  inputFile.close();
}

Motif::Motif() {
  // Default constructor implementation
}

Motif::Motif(std::string path) {
  std::string buffer;
  std::ifstream inputFile(path, std::ios::in);
  if (!inputFile.is_open()) {
    std::cout << "Motif file Open Failed " << std::endl;
    exit(1);
  } else {
    std::string line;
    std::getline(inputFile, line);
    std::istringstream iss(line);
    iss >> size >> type >> count;

    adjMatrix.resize(size, std::vector<ui>(size, 0));
    int x, y;
    while (std::getline(inputFile, line)) {
      std::istringstream iss(line);
      iss >> x >> y;
      adjMatrix[x][y] = 1;
      adjMatrix[y][x] = 1;
    }
  }

  inputFile.close();
  if (debug) {
    std::cout << "size=" << size << ", type=" << type << ", count=" << count
              << std::endl;
    cout << "adjMatrix " << endl;
    for (ui i = 0; i < size; i++) {
      for (ui j = 0; j < size; j++) {
        cout << adjMatrix[i][j] << " ";
      }
      cout << endl;
    }
  }
}
