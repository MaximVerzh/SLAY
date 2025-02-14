#ifndef CSR_MATRIX_H
#define CSR_MATRIX_H

#include <vector>
#include <map>
#include <stdexcept>

class CSRMatrix {
private:
    std::vector<double> values;
    std::vector<size_t> col_indices;
    std::vector<size_t> row_ptr;
    size_t rows, cols;

public:
    CSRMatrix(size_t rows, size_t cols, const std::map<std::pair<size_t, size_t>, double>& dok)
        : rows(rows), cols(cols), row_ptr(rows + 1, 0) {
        
        for (const auto& [pos, val] : dok) {
            if (val != 0) {
                values.push_back(val);
                col_indices.push_back(pos.second);
                row_ptr[pos.first + 1]++;
            }
        }

        for (size_t i = 1; i <= rows; i++) {
            row_ptr[i] += row_ptr[i - 1];
        }}

    double get(size_t row, size_t col) const {
        if (row >= rows || col >= cols) throw std::out_of_range("Index out of bounds");
        for (size_t i = row_ptr[row]; i < row_ptr[row + 1]; i++) {
            if (col_indices[i] == col) return values[i];
        }

        
        return 0.0;
    }

    size_t row_count() const {return rows;}
    size_t col_count() const { return cols;}
};

#endif
