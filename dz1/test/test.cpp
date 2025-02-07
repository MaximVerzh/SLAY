#include <gtest/gtest.h>
#include "solver.h"


TEST(TridiagonalSolverTest, SimpleCase) {
    std::vector<double> a = {0, -1, -1};
    std::vector<double> b = {2, 2, 2};
    std::vector<double> c = {-1, -1, 0};
    std::vector<double> d = {1, 2, 3};
    
    std::vector<double> expected = {2.5, 4, 3.5};
    auto result = TridiagonalSolver::solve(a, b, c, d);
    
    for (size_t i = 0; i < expected.size(); i++) {
        ASSERT_NEAR(result[i], expected[i], 1e-6);
    }
}

TEST(TridiagonalSolverTest, SimpleCase2) {
    std::vector<double> a = {0, -1, -1};
    std::vector<double> b = {3, 3, -2};
    std::vector<double> c = {-1, -1, 0};
    std::vector<double> d = {2, 3, 1};
    
    std::vector<double> expected = {1, 1, -1};
    auto result = TridiagonalSolver::solve(a, b, c, d);
    
    for (size_t i = 0; i < expected.size(); i++) {
        ASSERT_NEAR(result[i], expected[i], 1e-6);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
} 
