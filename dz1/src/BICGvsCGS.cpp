#include <iostream>
#include <fstream>
#include <vector>
#include "csr_matrix.h"
#include "bicg.h"
#include "cgs.h"
#include <cmath>
#include <random>


CSRMatrix generateRandomNonsymmetricMatrix(size_t size, double density = 0.1) {
    std::map<std::pair<size_t, size_t>, double> dok;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> val_dist(-1.0, 1.0);
    std::bernoulli_distribution bernoulli(density);

    for (size_t i = 0; i < size; ++i) {
        for (size_t j = 0; j < size; ++j) {
            if (bernoulli(gen)) {
                dok[{i, j}] = val_dist(gen);
            }
        }
    }
    
   
    for (size_t i = 0; i < size; ++i) {
        dok[{i, i}] = size * 2.0; 
    }

    return CSRMatrix(dok, size, size);
}


std::vector<double> generateRandomVector(size_t size) {
    std::vector<double> vec(size);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-10.0, 10.0);

    for (auto& elem : vec) {
        elem = dist(gen);
    }
    return vec;
}


std::pair<std::vector<double>, std::vector<double>> BiCGWithResidualHistory(
    const CSRMatrix& A, const std::vector<double>& b, size_t max_iter, double tolerance) {
    
    std::vector<double> x(b.size(), 0.0);
    std::vector<double> residuals;
    
    std::vector<double> r = b - A * x;
    std::vector<double> r_tilde = r;
    std::vector<double> p = r;
    std::vector<double> p_tilde = r_tilde;
    double rho = dotProduct(r_tilde, r);
    residuals.push_back(sqrt(dotProduct(r, r)));

    for (size_t iter = 0; iter < max_iter; ++iter) {
        std::vector<double> Ap = A * p;
        std::vector<double> ATp_tilde = A.transpose() * p_tilde;
        
        double alpha = rho / dotProduct(p_tilde, Ap);
        
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
            r_tilde[i] -= alpha * ATp_tilde[i];
        }

        double rho_new = dotProduct(r_tilde, r);
        residuals.push_back(sqrt(dotProduct(r, r)));

        if (residuals.back() < tolerance) break;

        double beta = rho_new / rho;
        for (size_t i = 0; i < p.size(); ++i) {
            p[i] = r[i] + beta * p[i];
            p_tilde[i] = r_tilde[i] + beta * p_tilde[i];
        }
        rho = rho_new;
    }

    return {x, residuals};
}

std::pair<std::vector<double>, std::vector<double>> CGSWithResidualHistory(
    const CSRMatrix& A, const std::vector<double>& b, size_t max_iter, double tolerance) {
    
    std::vector<double> x(b.size(), 0.0);
    std::vector<double> residuals;
    
    std::vector<double> r = b - A * x;
    std::vector<double> r_tilde = r;
    std::vector<double> u = r;
    std::vector<double> p = r;
    double rho = dotProduct(r_tilde, r);
    residuals.push_back(sqrt(dotProduct(r, r)));

    for (size_t iter = 0; iter < max_iter; ++iter) {
        std::vector<double> Ap = A * p;
        double alpha = rho / dotProduct(r_tilde, Ap);
        
        std::vector<double> q = u - alpha * Ap;
        std::vector<double> u_new = u + q;
        
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] += alpha * u_new[i];
        }

        std::vector<double> r_new = r - alpha * A * u_new;
        residuals.push_back(sqrt(dotProduct(r_new, r_new)));

        if (residuals.back() < tolerance) break;

        double rho_new = dotProduct(r_tilde, r_new);
        double beta = rho_new / rho;
        u = r_new + beta * q;
        p = u + beta * (q + beta * p);
        r = r_new;
        rho = rho_new;
    }

    return {x, residuals};
}


int main() {
    const size_t matrix_size = 250;
    const size_t max_iter = 1000;
    const double tolerance = 1e-8;
    
    
    CSRMatrix A = generateRandomNonsymmetricMatrix(matrix_size);
    std::vector<double> exact_solution = generateRandomVector(matrix_size);
    std::vector<double> b = A * exact_solution;
    
  
    auto [x_bicg, residuals_bicg] = BiCGWithResidualHistory(A, b, max_iter, tolerance);
    
   
    auto [x_cgs, residuals_cgs] = CGSWithResidualHistory(A, b, max_iter, tolerance);
    
 
    std::ofstream bicg_file("bicg_residuals.txt");
    std::ofstream cgs_file("cgs_residuals.txt");
    
    bicg_file << "Iteration\tResidual\n";
    for (size_t i = 0; i < residuals_bicg.size(); ++i) {
        bicg_file << i << "\t" << residuals_bicg[i] << "\n";
    }
    
    cgs_file << "Iteration\tResidual\n";
    for (size_t i = 0; i < residuals_cgs.size(); ++i) {
        cgs_file << i << "\t" << residuals_cgs[i] << "\n";
    }
    

    
    return 0;
}