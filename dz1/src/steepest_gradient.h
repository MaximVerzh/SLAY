#ifndef STEEPEST_GRADIENT_H
#define STEEPEST_GRADIENT_H

#include "csr_matrix.h"
#include <vector>
#include <cmath>
#include <iostream>

class SteepestGradient {
public:
    static std::vector<double> solve(
        const CSRMatrix& A,
        const std::vector<double>& b,
        double tolerance = 1e-6,
        size_t max_iterations = 1000
    ) {
        if (b.size() != A.getRows()) {
            throw std::invalid_argument("Размер вектора b должен совпадать с количеством строк матрицы A.");
        }

        size_t n = A.getRows();
        std::vector<double> x(n, 0.0); 
        std::vector<double> r = b;     

        for (size_t iter = 0; iter < max_iterations; iter++) {
            


            std::vector<double> Ar = A * r;

            
            double alpha_numerator = 0.0;
            double alpha_denominator = 0.0;
            for (size_t i = 0; i < n; i++) {
                alpha_numerator += r[i] * r[i];
                alpha_denominator += r[i] * Ar[i];
            }

            if (std::abs(alpha_denominator) < 1e-10) {
                std::cout << "Деление на ноль. Возможно, решение найдено.\n";
                return x;
            }

            double alpha = alpha_numerator / alpha_denominator;

           
            for (size_t i = 0; i < n; i++) {
                x[i] += alpha * r[i];
            }

           
            std::vector<double> Ax = A * x;
            for (size_t i = 0; i < n; i++) {
                r[i] = b[i] - Ax[i];
            }

            
            double residual_norm = 0.0;
            for (double val : r) {
                residual_norm += val * val;
            }
            residual_norm = std::sqrt(residual_norm);

            if (residual_norm < tolerance) {
                std::cout << "Наискорейший градиентный спуск: Сходимость достигнута за " 
                          << iter + 1 << " итераций.\n";
                return x;
            }
        }

        std::cout << "Наискорейший градиентный спуск: Достигнуто максимальное число итераций.\n";
        return x;
    }
};

#endif