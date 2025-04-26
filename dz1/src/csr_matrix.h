#ifndef CSRMATRIX_H
#define CSRMATRIX_H


#include <vector>
#include <map>
#include <stdexcept>
#include <algorithm>
#include <iostream> 
#include <cmath>    
#include <numeric>
class CSRMatrix {
private:
    std::vector<double> values;        
    std::vector<size_t> col_indices;   
    std::vector<size_t> row_ptr;       
    size_t rows, cols;                 

public:
    
    CSRMatrix(const std::map<std::pair<size_t, size_t>, double>& dok, size_t rows, size_t cols)
        : rows(rows), cols(cols) {
        if (dok.empty()) {
            throw std::invalid_argument("DOK is empty");
        }

        row_ptr.resize(rows + 1, 0);
        for (const auto& entry : dok) {
            size_t i = entry.first.first;
            size_t j = entry.first.second;
            double value = entry.second;

            if (i >= rows || j >= cols) {
                throw std::out_of_range("Index out of range in DOK");
            }

            values.push_back(value);
            col_indices.push_back(j);
            row_ptr[i + 1]++;
        }

        for (size_t i = 1; i <= rows; i++) {
            row_ptr[i] += row_ptr[i - 1];
        }
    }

  
    double operator()(size_t i, size_t j) const {
        if (i >= rows || j >= cols) {
            throw std::out_of_range("Index out of range");
        }

        for (size_t k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
            if (col_indices[k] == j) {
                return values[k];
            }
        }
        return 0.0;
    }

    
    std::vector<double> operator*(const std::vector<double>& vec) const {
        if (vec.size() != cols) {
            throw std::invalid_argument("Vector size must match matrix columns");
        }
        std::vector<double> result(rows, 0.0);
        for (size_t i = 0; i < rows; i++) {
            for (size_t k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
                result[i] += values[k] * vec[col_indices[k]];
            }
        }
        return result;
    }

  
    std::vector<double> jacobi(const std::vector<double>& b, double tolerance = 1e-6, size_t max_iterations = 1000) const {
        if (b.size() != rows) {
            throw std::invalid_argument("Vector size must match matrix rows");
        }

  
        std::vector<double> x(rows, 0.0);

     
        for (size_t iter = 0; iter < max_iterations; iter++) {
            std::vector<double> x_new(rows, 0.0); 

            for (size_t i = 0; i < rows; i++) {
                double sum = 0.0;

               
                for (size_t k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
                    size_t j = col_indices[k];
                    if (j != i) { 
                        sum += values[k] * x[j];
                    }
                }

          
                double diag = (*this)(i, i);
                if (diag == 0.0) {
                    throw std::runtime_error("Zero diagonal element detected");
                }
                x_new[i] = (b[i] - sum) / diag;
            }


            double diff = 0.0;
            for (size_t i = 0; i < rows; i++) {
                diff += std::abs(x_new[i] - x[i]);
            }

            if (diff < tolerance) {
                std::cout << "Метод Якоби: Итерации сошлись за " << iter + 1 << " шагов.\n";
                return x_new;
            }

         
            x = x_new;
        }

        std::cout << "Метод Якоби: Достигнуто максимальное число итераций.\n";
        return x;
    }


    std::vector<double> gaussSeidel(const std::vector<double>& b, double tolerance = 1e-6, size_t max_iterations = 1000) const {
        if (b.size() != rows) {
            throw std::invalid_argument("Vector size must match matrix rows");
        }
        std::vector<double> x(rows, 0.0);

        for (size_t iter = 0; iter < max_iterations; iter++) {
            std::vector<double> x_new = x; 

            for (size_t i = 0; i < rows; i++) {
                double sum = 0.0;

             
                for (size_t k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
                    size_t j = col_indices[k];
                    if (j != i) { 
                        sum += values[k] * x_new[j]; 
                    }
                }

               
                double diag = (*this)(i, i);
                if (diag == 0.0) {
                    throw std::runtime_error("Zero diagonal element detected");
                }
                x_new[i] = (b[i] - sum) / diag;
            }


            double diff = 0.0;
            for (size_t i = 0; i < rows; i++) {
                diff += std::abs(x_new[i] - x[i]);
            }

            if (diff < tolerance) {
                std::cout << "Метод Гаусса-Зейделя: Итерации сошлись за " << iter + 1 << " шагов.\n";
                return x_new;
            }

          
            x = x_new;
        }

        std::cout << "Метод Гаусса-Зейделя: Достигнуто максимальное число итераций.\n";
        return x;
    }

   
    
    std::vector<double> simpleIteration(const std::vector<double>& b, double tau = 0.1, double tolerance = 1e-6, size_t max_iterations = 1000) const {
        if (b.size() != rows) {
            throw std::invalid_argument("Vector size must match matrix rows");
        }

      
        std::vector<double> x(rows, 0.0);

    
        for (size_t iter = 0; iter < max_iterations; iter++) {
            std::vector<double> Ax = (*this) * x; 
            std::vector<double> x_new(rows);

  
            for (size_t i = 0; i < rows; i++) {
                x_new[i] = x[i] - tau * (Ax[i] - b[i]);
            }
            double diff = 0.0;
            for (size_t i = 0; i < rows; i++) {
                diff += std::abs(x_new[i] - x[i]);
            }

            if (diff < tolerance) {
                std::cout << "Метод простой итерации: Итерации сошлись за " << iter + 1 << " шагов.\n";
                return x_new;
            }

            x = x_new;
        }

        std::cout << "Метод простой итерации: Достигнуто максимальное число итераций.\n";
        return x;
    }
    double powerMethod(double tolerance = 1e-6, size_t max_iterations = 1000) const {
        
        std::vector<double> r(rows, 1.0);
        double lambda_old = 0.0; 
        double lambda_new = 0.0; 
    
        for (size_t iter = 0; iter < max_iterations; iter++) {
          
            std::vector<double> r_new = (*this) * r;
            double norm = 0.0;
            for (double val : r_new) {
                norm += val * val;
            }
            norm = std::sqrt(norm);
    
    
            for (size_t i = 0; i < rows; i++) {
                r_new[i] /= norm;
            }
    
           

            std::vector<double> Ar_new = (*this) * r_new; 
            double numerator = 0.0;
            double denominator = 0.0;
            for (size_t i = 0; i < rows; i++) {
                numerator += r_new[i] * Ar_new[i]; 
                denominator += r_new[i] * r_new[i];
            }
            lambda_new = numerator / denominator;
    
          
            if (std::abs(lambda_new - lambda_old) < tolerance) {
                std::cout << "Степенной метод: Сходимость достигнута за " << iter + 1 << " итераций.\n";
                return lambda_new;
            }
    
           
            r = r_new;
            lambda_old = lambda_new;
        }
    
        std::cout << "Степенной метод: Достигнуто максимальное число итераций.\n";
        return lambda_new;
    }
    std::vector<double> chebyshevAcceleration(const std::vector<double>& b, double lambda_min, double lambda_max, size_t n, double tolerance = 1e-6, size_t max_iterations = 1000) const {
        if (b.size() != rows) {
            throw std::invalid_argument("Vector size must match matrix rows");
        }

        std::vector<double> x(rows, 0.0);
        std::vector<double> tau_k(n);

       
        for (size_t k = 0; k < n; k++) {
            double t_k = cos(M_PI * (2 * k + 1) / (2 * n));
            tau_k[k] = 2 / (lambda_max + lambda_min + (lambda_max - lambda_min) * t_k);
        }

        for (size_t iter = 0; iter < max_iterations; iter++) {
            std::vector<double> x_new = x;

            for (size_t k = 0; k < n; k++) {
                std::vector<double> Ax = (*this) * x_new;
                for (size_t i = 0; i < rows; i++) {
                    x_new[i] = x_new[i] - tau_k[k] * (Ax[i] - b[i]);
                }
            }

            double diff = 0.0;
            for (size_t i = 0; i < rows; i++) {
                diff += std::abs(x_new[i] - x[i]);
            }

            if (diff < tolerance) {
                std::cout << "Чебышевское ускорение: Итерации сошлись за " << iter + 1 << " шагов.\n";
                return x_new;
            }

            x = x_new;
        }

        std::cout << "Чебышевское ускорение: Достигнуто максимальное число итераций.\n";
        return x;
    }
    CSRMatrix transpose() const {
        std::map<std::pair<size_t, size_t>, double> dok;
        
        for (size_t i = 0; i < rows; i++) {
            for (size_t k = row_ptr[i]; k < row_ptr[i + 1]; k++) {
                size_t j = col_indices[k];
                dok[{j, i}] = values[k]; 
            }
        }
        
        return CSRMatrix(dok, cols, rows); 
    }
    size_t solveWithIterCount(
        const std::vector<double>& b, 
        double tolerance, 
        size_t max_iter,
        const std::string& method
    ) const {
        std::vector<double> x(rows, 0.0);
        for (size_t iter = 0; iter < max_iter; iter++) {
           
        }
        return max_iter;
    }
    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }
    const std::vector<double>& getValues() const { return values; }
    const std::vector<size_t>& getColIndices() const { return col_indices; }
    const std::vector<size_t>& getRowPtr() const { return row_ptr; }
};

inline std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must have same size for subtraction");
    }
    
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

inline std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must have same size for addition");
    }
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

inline std::vector<double> operator*(double scalar, const std::vector<double>& vec) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = scalar * vec[i];
    }
    return result;
}
inline CSRMatrix operator*(double scalar, const CSRMatrix& mat) {
    std::map<std::pair<size_t, size_t>, double> dok;
    
    for (size_t i = 0; i < mat.getRows(); i++) {
        for (size_t k = mat.getRowPtr()[i]; k < mat.getRowPtr()[i+1]; k++) {
            size_t j = mat.getColIndices()[k];
            dok[{i, j}] = mat.getValues()[k] * scalar;
        }
    }
    
    return CSRMatrix(dok, mat.getRows(), mat.getCols());
}
inline double dotProduct(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must have same size for dot product");
    }
    double result = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        result += a[i] * b[i];
    }
    return result;
}
#endif               