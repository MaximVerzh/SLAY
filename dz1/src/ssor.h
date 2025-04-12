#ifndef SSOR_H
#define SSOR_H

#include "csr_matrix.h"
#include <vector>
#include <cmath>
#include <iostream>

class SSOR {
public:
    static std::vector<double> solve(
        const CSRMatrix& A,
        const std::vector<double>& b,
        double omega = 1.0,  
        double tolerance = 1e-6,
        size_t max_iterations = 1000
    ) {
        if (b.size() != A.getRows()) {
            throw std::invalid_argument("Размер вектора b должен совпадать с количеством строк матрицы A.");
        }

        size_t n = A.getRows();
        std::vector<double> x(n, 0.0);  

        for (size_t iter = 0; iter < max_iterations; iter++) {
            std::vector<double> x_prev = x;

          
            for (size_t i = 0; i < n; i++) {
                double sum = 0.0;
                double diag = 0.0;

                for (size_t k = A.getRowPtr()[i]; k < A.getRowPtr()[i + 1]; k++) {
                    size_t j = A.getColIndices()[k];
                    if (j != i) {
                        sum += A.getValues()[k] * x[j];
                    } else {
                        diag = A.getValues()[k];
                    }
                }

                if (diag == 0.0) {
                    throw std::runtime_error("Нулевой диагональный элемент в матрице.");
                }

                x[i] = (1 - omega) * x[i] + (omega / diag) * (b[i] - sum);
            }

           
            for (int i = n - 1; i >= 0; i--) {
                double sum = 0.0;
                double diag = 0.0;

                for (size_t k = A.getRowPtr()[i]; k < A.getRowPtr()[i + 1]; k++) {
                    size_t j = A.getColIndices()[k];
                    if (j != i) {
                        sum += A.getValues()[k] * x[j];
                    } else {
                        diag = A.getValues()[k];
                    }
                }

                if (diag == 0.0) {
                    throw std::runtime_error("Нулевой диагональный элемент в матрице.");
                }

                x[i] = (1 - omega) * x[i] + (omega / diag) * (b[i] - sum);
            }

           
            double diff = 0.0;
            for (size_t i = 0; i < n; i++) {
                diff += std::abs(x[i] - x_prev[i]);
            }

            if (diff < tolerance) {
                std::cout << "SSOR: Сходимость достигнута за " << iter + 1 << " итераций.\n";
                return x;
            }
        }

        std::cout << "SSOR: Достигнуто максимальное число итераций.\n";
        return x;
    }
};

#endif