#ifndef CSRMATRIX_H
#define CSRMATRIX_H


#include <vector>
#include <map>
#include <stdexcept>
#include <algorithm>
#include <iostream> 
#include <cmath>    

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


    size_t getRows() const { return rows; }
    size_t getCols() const { return cols; }
    const std::vector<double>& getValues() const { return values; }
    const std::vector<size_t>& getColIndices() const { return col_indices; }
    const std::vector<size_t>& getRowPtr() const { return row_ptr; }
};
#endif               