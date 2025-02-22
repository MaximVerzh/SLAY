#ifndef VECTOR_OPS_H
#define VECTOR_OPS_H

#include <vector>
#include <stdexcept>


inline std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] + b[i];
    }
    return result;
}


inline double operator*(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }
    double result = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        result += a[i] * b[i];
    }
    return result;
}


inline std::vector<double> operator*(const std::vector<double>& a, double scalar) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] * scalar;
    }
    return result;
}


inline std::vector<double> operator*(double scalar, const std::vector<double>& a) {
    return a * scalar;
}

#endif 