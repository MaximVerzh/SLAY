#ifndef CHEBYSHEV_SYMMETRIC_H
#define CHEBYSHEV_SYMMETRIC_H

#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <iostream>

template <typename Vector>
class ChebyshevSymmetricSolver {
private:
    double rho_;
    size_t max_iter_;
    double tolerance_;
    static double norm_diff(const Vector& a, const Vector& b) {
        double s = 0.0;
        for (size_t i = 0; i < a.size(); ++i) {
            s += std::abs(a[i] - b[i]);
        }
        return s;
    }

public:

    ChebyshevSymmetricSolver(double rho, size_t max_iter = 1000, double tolerance = 1e-6)
        : rho_(rho), max_iter_(max_iter), tolerance_(tolerance) {
        if (rho <= 0.0 || rho >= 1.0) {
            throw std::invalid_argument("Спектральный радиус rho должен быть в (0, 1)");
        }
    }

  
    Vector solve(const Vector& x0, std::function<Vector(const Vector&)> step_func) const {
        size_t n = x0.size();
        Vector y_prev = x0;
        Vector y_curr = step_func(y_prev);
        Vector y_next(n);

        double omega = 2.0 / (2.0 - rho_ * rho_);
        double diff = 0.0;

        for (size_t iter = 1; iter < max_iter_; ++iter) {
         
            Vector tmp = step_func(y_curr);
            for (size_t i = 0; i < n; ++i) {
                y_next[i] = omega * (tmp[i] - y_prev[i]) + y_prev[i];
            }

            diff = norm_diff(y_next, y_curr);
            if (diff < tolerance_) {
                std::cout << "Чебышёвское ускорение: Сходимость достигнута за " 
                          << iter + 1 << " итераций.\n";
                return y_next;
            }

         
            omega = 1.0 / (1.0 - (rho_ * rho_ * omega) / 4.0);
            y_prev = y_curr;
            y_curr = y_next;
        }

        std::cout << "Чебышёвское ускорение: Достигнуто максимальное число итераций.\n";
        return y_next;
    }

    
    void set_rho(double rho) {
        if (rho <= 0.0 || rho >= 1.0) {
            throw std::invalid_argument("Спектральный радиус rho должен быть в (0, 1)");
        }
        rho_ = rho;
    }

    void set_max_iter(size_t max_iter) {
        max_iter_ = max_iter;
    }

    void set_tolerance(double tolerance) {
        tolerance_ = tolerance;
    }
};

#endif
