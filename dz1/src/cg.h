#ifndef CG_H
#define CG_H

#include "csr_matrix.h"
#include <vector>
#include <cmath>

class ConjugateGradient {
public:
    static std::vector<double> solve(const CSRMatrix& A, 
                                   const std::vector<double>& b,
                                   size_t max_iter = 1000,
                                   double tolerance = 1e-6,
                                   const std::vector<double>& x0 = {}) {
        
        if (b.size() != A.getRows()) {
            throw std::invalid_argument("Matrix and vector dimensions mismatch");
        }

        std::vector<double> x = x0.empty() ? std::vector<double>(b.size(), 0.0) : x0;
        std::vector<double> r = A * x - b;
        std::vector<double> d = r;
        double delta_new = dotProduct(r, r);
        const double delta0 = delta_new;

        for (size_t iter = 0; iter < max_iter; ++iter) {
            std::vector<double> Ad = A * d;
            double alpha = delta_new / dotProduct(d, Ad);
            
            for (size_t i = 0; i < x.size(); ++i) {
                x[i] -= alpha * d[i];
                r[i] -= alpha * Ad[i];
            }

            double delta_old = delta_new;
            delta_new = dotProduct(r, r);
            
            if (std::sqrt(delta_new) < tolerance * std::sqrt(delta0)) {
                std::cout << "CG converged in " << iter + 1 << " iterations\n";
                return x;
            }

            double beta = delta_new / delta_old;
            for (size_t i = 0; i < d.size(); ++i) {
                d[i] = r[i] + beta * d[i];
            }
        }

        std::cout << "CG reached maximum iterations\n";
        return x;
    }

private:
    static double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
        double result = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            result += a[i] * b[i];
        }
        return result;
    }
};

#endif