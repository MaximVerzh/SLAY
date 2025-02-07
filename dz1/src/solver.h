#ifndef SOLVER_H
#define SOLVER_H

#include <vector>

class TridiagonalSolver {
public:
    static std::vector<double> solve(const std::vector<double>& a, const std::vector<double>& b,
                                     const std::vector<double>& c, const std::vector<double>& d) {
        int n = b.size();
        std::vector<double> p(n, 0), q(n, 0), x(n, 0);
        
        p[0] = c[0] / b[0];
        q[0] = d[0] / b[0];
        
        for (int i = 1; i < n; i++) {
            double denom = b[i] - a[i] * p[i - 1];
            p[i] = c[i] / denom;
            q[i] = (d[i] - a[i] * q[i - 1]) / denom;
        }
        
        x[n - 1] = q[n - 1];
        for (int i = n - 2; i >= 0; i--) {
            x[i] = q[i] - p[i] * x[i + 1];
        }
        
        return x;
    }
};

#endif 
