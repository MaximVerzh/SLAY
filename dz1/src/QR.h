#ifndef QR_DECOMPOSITION_H
#define QR_DECOMPOSITION_H

#include "matrix.h"
#include <vector>
#include <cmath>
#include <stdexcept>

void householderQR(const Matrix& A, Matrix& Q, Matrix& R) {
    size_t m = A.getRows();
    size_t n = A.getCols();

    if (m < n) {
        throw std::invalid_argument("Number of rows must be greater than or equal to number of columns");
    }

    Q = Matrix(m, m);
    R = A;

    for (size_t i = 0; i < m; i++) {
        Q(i, i) = 1.0;
    }

    for (size_t k = 0; k < n; k++) {
        double norm = 0.0;
        for (size_t i = k; i < m; i++) {
            norm += R(i, k) * R(i, k);
        }
        norm = std::sqrt(norm);

        double alpha = -std::copysign(norm, R(k, k));
        double u1 = R(k, k) - alpha;
        double norm_u = std::sqrt(u1 * u1 + norm * norm - R(k, k) * R(k, k));

        std::vector<double> v(m, 0.0);
        v[k] = u1 / norm_u;
        for (size_t i = k + 1; i < m; i++) {
            v[i] = R(i, k) / norm_u;
        }

        for (size_t j = k; j < n; j++) {
            double dot = 0.0;
            for (size_t i = k; i < m; i++) {
                dot += v[i] * R(i, j);
            }
            for (size_t i = k; i < m; i++) {
                R(i, j) -= 2.0 * v[i] * dot;
            }
        }

        for (size_t i = k + 1; i < m; i++) {
            R(i, k) = 0.0;
        }

        for (size_t j = 0; j < m; j++) {
            double dot = 0.0;
            for (size_t i = k; i < m; i++) {
                dot += v[i] * Q(j, i);
            }
            for (size_t i = k; i < m; i++) {
                Q(j, i) -= 2.0 * v[i] * dot;
            }
        }
    }
}

std::vector<double> solveQR(const Matrix& A, const std::vector<double>& b) {
    size_t m = A.getRows();
    size_t n = A.getCols();

    if (m != b.size()) {
        throw std::invalid_argument("Matrix and vector dimensions must match");
    }

    Matrix Q(m, m), R(m, n);
    householderQR(A, Q, R);

    std::vector<double> Qt_b(m, 0.0);
    for (size_t i = 0; i < m; i++) {
        for (size_t j = 0; j < m; j++) {
            Qt_b[i] += Q(j, i) * b[j];
        }
    }

    std::vector<double> x(n, 0.0);
    for (int i = n - 1; i >= 0; i--) {
        x[i] = Qt_b[i];
        for (size_t j = i + 1; j < n; j++) {
            x[i] -= R(i, j) * x[j];
        }
        x[i] /= R(i, i);
    }

    return x;
}

#endif // QR_DECOMPOSITION_H