#ifndef BICG_H
#define BICG_H

#include "csr_matrix.h"
#include <vector>
#include <cmath>

class BiConjugateGradient {
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
        std::vector<double> r_tilde = r;
        std::vector<double> p = r;
        std::vector<double> p_tilde = r_tilde;
        double rho = dotProduct(r_tilde, r);

        for (size_t iter = 0; iter < max_iter; ++iter) {
            std::vector<double> Ap = A * p;
            std::vector<double> ATp_tilde = A.transpose() * p_tilde;
            
            double alpha = rho / dotProduct(p_tilde, Ap);
            
            for (size_t i = 0; i < x.size(); ++i) {
                x[i] -= alpha * p[i];
                r[i] -= alpha * Ap[i];
                r_tilde[i] -= alpha * ATp_tilde[i];
            }

            double rho_new = dotProduct(r_tilde, r);
            if (std::sqrt(dotProduct(r, r)) < tolerance) {
                std::cout << "BiCG converged in " << iter + 1 << " iterations\n";
                return x;
            }

            double beta = rho_new / rho;
            for (size_t i = 0; i < p.size(); ++i) {
                p[i] = r[i] + beta * p[i];
                p_tilde[i] = r_tilde[i] + beta * p_tilde[i];
            }
            rho = rho_new;
        }

        std::cout << "BiCG reached maximum iterations\n";
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
