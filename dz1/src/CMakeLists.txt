cmake_minimum_required(VERSION 3.10)
project(TridiagonalSolver)

add_library(solver INTERFACE)
target_include_directories(solver INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})