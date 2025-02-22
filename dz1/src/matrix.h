#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <stdexcept>

class Matrix {
    private:
        std::vector<double> data;
        size_t rows, cols;
    
    public:
    
        Matrix(size_t rows, size_t cols) : rows(rows), cols(cols), data(rows * cols, 0.0) {}
    
        
        Matrix(size_t rows, size_t cols, const std::vector<double>& values)
            : rows(rows), cols(cols), data(values) {
            if (values.size() != rows * cols) {
                throw std::invalid_argument("Number of values must match matrix dimensions");
            }
        }
    
      
        double& operator()(size_t i, size_t j) {
            if (i >= rows || j >= cols) {
                throw std::out_of_range("Index out of range");
            }
            return data[i * cols + j];
        }
    
        const double& operator()(size_t i, size_t j) const {
            if (i >= rows || j >= cols) {
                throw std::out_of_range("Index out of range");
            }
            return data[i * cols + j];
        }
    
        std::vector<double> operator*(const std::vector<double>& vec) const {
            if (vec.size() != cols) {
                throw std::invalid_argument("Vector size must match matrix columns");
            }
            std::vector<double> result(rows, 0.0);
            for (size_t i = 0; i < rows; i++) {
                for (size_t j = 0; j < cols; j++) {
                    result[i] += (*this)(i, j) * vec[j];
                }
            }
            return result;
        }
    
        
    
        size_t getRows() const { return rows; }
        size_t getCols() const { return cols; }
    };


#endif 