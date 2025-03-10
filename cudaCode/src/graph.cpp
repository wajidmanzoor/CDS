#include "../inc/graph.h"

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

        offset.resize(n + 1, 0);
        neighbors.resize(2 * m);
        degree.resize(n);
        int vertex, neigh;
        while (std::getline(inputFile, line)) {
            std::istringstream iss(line);
            iss >> vertex;
            while (iss >> neigh) {
                if (vertex == neigh) continue;
                neighbors[offset[vertex] + offset[vertex + 1]] = neigh;
                offset[vertex + 1]++;
            }
            degree[vertex] = offset[vertex + 1];
        }
    }

    inputFile.close();
    std::cout << "n =" << n << ", m=" << m << std::endl;
}

Graph::Graph() {
    // Default constructor implementation
}

void Graph::printGraph() {
    std::cout << "Print Hello" << std::endl;
}

void Graph::getListingOrder(std::vector<ui>& arr) {
    std::vector<ui> sortedbyCore;
    coreDecompose(sortedbyCore);

    for (size_t i = 0; i < n; ++i) {
        arr[corePeelSequence[i]] = i + 1;
    }
}

void Graph::coreDecompose(std::vector<ui>& arr) {
    core.resize(n);
    int maxDegree = *std::max_element(degree.begin(), degree.end());

    std::vector<ui> bins(maxDegree + 1, 0);

    for (ui deg : degree) {
        bins[deg]++;
    }

    std::vector<int> bin_positions(maxDegree + 1, 0);
    std::partial_sum(bins.begin(), bins.end(), bin_positions.begin());

    std::vector<ui> position(n + 1);
    std::vector<ui> sortedVertex(n + 1);

    for (ui v = 0; v < n; v++) {
        position[v] = bins[degree[v]];
        sortedVertex[position[v]] = v;
        bins[degree[v]]++;
    }

    for (int i = maxDegree; i >= 1; i--) {
        bins[i] = bins[i - 1];
    }

    bins[0] = 1;

    for (int i = 0; i < n; i++) {
        ui v = sortedVertex[i];
        for (int j = offset[v]; j < offset[v + 1]; j++) {
            ui u = neighbors[j];
            if (degree[u] > degree[v]) {
                ui du = degree[u];
                ui pu = position[u];
                ui pw = bins[du];
                ui w = sortedVertex[pw];
                if (u != w) {
                    position[u] = pw;
                    sortedVertex[pu] = w;
                    position[w] = pu;
                    sortedVertex[pw] = u;
                }

                bins[du]++;
                degree[u]--;
            }
        }

        arr[n - i - 1] = v;
    }
}
