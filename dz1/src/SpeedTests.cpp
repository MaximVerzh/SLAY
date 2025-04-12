#include "csr_matrix.h"
#include "sor.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>


std::map<std::pair<size_t, size_t>, double> generate_spd_matrix(size_t n) {
    std::map<std::pair<size_t, size_t>, double> dok;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(0.1, 1.0);

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j <= i; ++j) {
            double val = dist(gen);
            dok[{i, j}] = val;
            if (i != j) dok[{j, i}] = val;
        }
        dok[{i, i}] += n; 
    }
    return dok;
}


void write_error_history(const std::vector<double>& errors, const std::string& filename) {
    std::ofstream out(filename);
    out << "iteration,error\n";
    for (size_t i = 0; i < errors.size(); ++i) {
        out << i+1 << "," << errors[i] << "\n";
    }
    out.close();
}


std::vector<double> jacobiWithHistory(const CSRMatrix& A, const std::vector<double>& b, 
                                     double tolerance, size_t max_iter, 
                                     const std::vector<double>& exact_solution) {
    std::vector<double> errors;
    std::vector<double> x(A.getRows(), 0.0);
    
    for (size_t iter = 0; iter < max_iter; ++iter) {
        std::vector<double> x_new = x;
        
        for (size_t i = 0; i < A.getRows(); ++i) {
            double sum = 0.0;
            double diag = 0.0;
            
            for (size_t k = A.getRowPtr()[i]; k < A.getRowPtr()[i+1]; ++k) {
                size_t j = A.getColIndices()[k];
                if (j != i) sum += A.getValues()[k] * x[j];
                else diag = A.getValues()[k];
            }
            
            x_new[i] = (b[i] - sum) / diag;
        }
        
       
        double error = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            error += std::abs(x_new[i] - exact_solution[i]);
        }
        errors.push_back(error);
        
        x = x_new;
        
        if (error < tolerance) break;
    }
    
    return errors;
}



std::vector<double> gaussSeidelWithHistory(const CSRMatrix& A, const std::vector<double>& b, 
                                          double tolerance, size_t max_iter,
                                          const std::vector<double>& exact_solution) {
    std::vector<double> errors;
    std::vector<double> x(A.getRows(), 0.0);
    
    for (size_t iter = 0; iter < max_iter; ++iter) {
        std::vector<double> x_old = x;
        
        for (size_t i = 0; i < A.getRows(); ++i) {
            double sum = 0.0;
            double diag = 0.0;
            
            for (size_t k = A.getRowPtr()[i]; k < A.getRowPtr()[i+1]; ++k) {
                size_t j = A.getColIndices()[k];
                if (j != i) sum += A.getValues()[k] * x[j];
                else diag = A.getValues()[k];
            }
            
            x[i] = (b[i] - sum) / diag;
        }
        
     
        double error = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            error += std::abs(x[i] - exact_solution[i]);
        }
        errors.push_back(error);
        
        if (error < tolerance) break;
    }
    
    return errors;
}




std::vector<double> sorWithHistory(const CSRMatrix& A, const std::vector<double>& b, 
                                  double omega, double tolerance, size_t max_iter,
                                  const std::vector<double>& exact_solution) {
    std::vector<double> errors;
    std::vector<double> x(A.getRows(), 0.0);
    
    for (size_t iter = 0; iter < max_iter; ++iter) {
        std::vector<double> x_old = x;
        
        for (size_t i = 0; i < A.getRows(); ++i) {
            double sum = 0.0;
            double diag = 0.0;
            
            for (size_t k = A.getRowPtr()[i]; k < A.getRowPtr()[i+1]; ++k) {
                size_t j = A.getColIndices()[k];
                if (j != i) sum += A.getValues()[k] * x[j];
                else diag = A.getValues()[k];
            }
            
            x[i] = (1 - omega) * x_old[i] + (omega / diag) * (b[i] - sum);
        }
        
       
        double error = 0.0;
        for (size_t i = 0; i < x.size(); ++i) {
            error += std::abs(x[i] - exact_solution[i]);
        }
        errors.push_back(error);
        
        if (error < tolerance) break;
    }
    
    return errors;
}

int main() {
    const size_t n = 10;
    auto dok = generate_spd_matrix(n);
    CSRMatrix A(dok, n, n);



    std::vector<double> exact_solution(n);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> dist(-10.0, 10.0);
    for (auto& val : exact_solution) val = dist(gen);

    std::vector<double> b = A * exact_solution;
    const double tolerance = 1e-10;
    const size_t max_iter = 1000;
    const double omega = 1.02;


    auto errors_jacobi = jacobiWithHistory(A, b, tolerance, max_iter, exact_solution);
    auto errors_gs = gaussSeidelWithHistory(A, b, tolerance, max_iter, exact_solution);
    auto errors_sor = sorWithHistory(A, b, omega, tolerance, max_iter, exact_solution);

   
    write_error_history(errors_jacobi, "jacobi_errors.csv");
    write_error_history(errors_gs, "gs_errors.csv");
    write_error_history(errors_sor, "sor_errors.csv");

 

    return 0;
}