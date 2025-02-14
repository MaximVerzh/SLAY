#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <stdexcept>

class Matrix {
private:
    std::vector<std::vector<double>> data;
    size_t rows, cols;

public:


    Matrix(size_t rows, size_t cols, double value = 0.0) 
        : rows(rows), cols(cols), data(rows, std::vector<double>(cols, value)) {}

    double& get(size_t row, size_t col) {
        if (row >= rows || col >= cols) throw std::out_of_range("Index out of bounds");
        return data[row][col];
    }

    const double& get(size_t row, size_t col) const {
        if (row >= rows || col >= cols) throw std::out_of_range("Index out of bounds");
        return data[row][col];
}

    size_t row_count() const {return rows;}
    size_t col_count() const {return cols;}
};

#endif