#ifndef ELLIPTIC_MATRIX_GENERATOR_H
#define ELLIPTIC_MATRIX_GENERATOR_H

#include "csr_matrix.h"
#include <map>
#include <cmath>

class EllipticMatrixGenerator {
public:
    static CSRMatrix generatePoissonMatrix(size_t grid_size) {
        size_t N = grid_size;
        size_t total_nodes = N * N;
        std::map<std::pair<size_t, size_t>, double> dok;
        double h = 1.0 / (N - 1);

        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                size_t node = i * N + j;

               
                dok[{node, node}] = 4.0 / (h * h);

              
                if (j > 0)   dok[{node, node - 1}] = -1.0 / (h * h); // Left
                if (j < N-1) dok[{node, node + 1}] = -1.0 / (h * h); // Right
                if (i > 0)   dok[{node, node - N}] = -1.0 / (h * h); // Top
                if (i < N-1) dok[{node, node + N}] = -1.0 / (h * h); // Bottom
            }
        }

        return CSRMatrix(dok, total_nodes, total_nodes);
    }

    static std::vector<double> generateRHS(size_t grid_size) {
        size_t N = grid_size;
        std::vector<double> b(N * N);
        double h = 1.0 / (N - 1);

        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                double x = j * h;
                double y = i * h;
                b[i * N + j] = sin(M_PI * x) * sin(M_PI * y);
            }
        }

        return b;
    }
};

#endif