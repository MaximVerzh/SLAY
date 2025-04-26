#ifndef CGS_H
#define CGS_H

#include "csr_matrix.h"
#include <vector>
#include <cmath>

class ConjugateGradientSquared {
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
        std::vector<double> u = r;
        std::vector<double> p = r;
        double rho = dotProduct(r_tilde, r);

        for (size_t iter = 0; iter < max_iter; ++iter) {
            std::vector<double> Ap = A * p;
            double alpha = rho / dotProduct(r_tilde, Ap);
            
            std::vector<double> q = u - alpha * Ap;
            std::vector<double> u_new = u + q;
            
            for (size_t i = 0; i < x.size(); ++i) {
                x[i] += alpha * u_new[i];
            }

            std::vector<double> r_new = r - (A * (alpha * u_new));
            if (std::sqrt(dotProduct(r_new, r_new)) < tolerance) {
                std::cout << "CGS converged in " << iter + 1 << " iterations\n";
                return x;
            }

            double rho_new = dotProduct(r_tilde, r_new);
            double beta = rho_new / rho;
            u = r_new + beta * q;
            p = u + beta * (q + beta * p);
            r = r_new;
            rho = rho_new;
        }

        std::cout << "CGS reached maximum iterations\n";
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