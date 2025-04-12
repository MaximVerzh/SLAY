#ifndef SOR_H
#define SOR_H

#include "csr_matrix.h"
#include <vector>
#include <cmath>
#include <iostream>

class SOR {
public:
    static std::vector<double> solve(
        const CSRMatrix& A,
        const std::vector<double>& b,
        double omega = 1.5,  
        double tolerance = 1e-6,
        size_t max_iterations = 1000
    ) {
        if (omega <= 0 || omega >= 2) {
            throw std::invalid_argument("Параметр релаксации ω должен быть в интервале (0, 2)");
        }

        size_t n = A.getRows();
        std::vector<double> x(n, 0.0);

        for (size_t iter = 0; iter < max_iterations; iter++) {
            std::vector<double> x_prev = x;

            for (size_t i = 0; i < n; i++) {
                double sigma = 0.0;
                double diag = 0.0;

              
                for (size_t k = A.getRowPtr()[i]; k < A.getRowPtr()[i + 1]; k++) {
                    size_t j = A.getColIndices()[k];
                    if (j != i) {
                        sigma += A.getValues()[k] * x[j];
                    } else {
                        diag = A.getValues()[k];
                    }
                }

                if (diag == 0.0) {
                    throw std::runtime_error("Нулевой диагональный элемент в матрице");
                }

               
                x[i] = (1 - omega) * x_prev[i] + (omega / diag) * (b[i] - sigma);
            }

        


            double diff = 0.0;
            for (size_t i = 0; i < n; i++) {
                diff += std::abs(x[i] - x_prev[i]);
            }

            if (diff < tolerance) {
                std::cout << "SOR: Сходимость достигнута за " << iter + 1 << " итераций (ω=" << omega << ")\n";
                return x;
            }
        }

        std::cout << "SOR: Достигнуто максимальное число итераций\n";
        return x;
    }
    static size_t solveWithIterCount(
        const CSRMatrix& A,
        const std::vector<double>& b,
        double omega,
        double tolerance = 1e-6,
        size_t max_iterations = 1000
    ) {
        size_t n = A.getRows();
        std::vector<double> x(n, 0.0);
    
        for (size_t iter = 0; iter < max_iterations; iter++) {
            std::vector<double> x_prev = x;
    
            for (size_t i = 0; i < n; i++) {
                double sigma = 0.0;
                double diag = 0.0;
    
                for (size_t k = A.getRowPtr()[i]; k < A.getRowPtr()[i + 1]; k++) {
                    size_t j = A.getColIndices()[k];
                    if (j != i) {
                        sigma += A.getValues()[k] * x[j];
                    } else {
                        diag = A.getValues()[k];
                    }
                }
    
                x[i] = (1 - omega) * x_prev[i] + (omega / diag) * (b[i] - sigma);
            }
    
            double diff = 0.0;
            for (size_t i = 0; i < n; i++) {
                diff += std::abs(x[i] - x_prev[i]);
            }
    
            if (diff < tolerance) {
                return iter + 1;  // Возвращаем количество итераций
            }
        }
    
        return max_iterations;
    }
};

#endif