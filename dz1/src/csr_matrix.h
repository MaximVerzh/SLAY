#ifndef CSRMATRIX_H
#define CSRMATRIX_H

#include <vector>
#include <map>
#include <stdexcept>
#include <algorithm>

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
   
    size_t getRows() const { return rows; }

    
    size_t getCols() const { return cols; }

   
    const std::vector<double>& getValues() const { return values; }

   
    const std::vector<size_t>& getColIndices() const { return col_indices; }

   
    const std::vector<size_t>& getRowPtr() const { return row_ptr; }
};

#endif 