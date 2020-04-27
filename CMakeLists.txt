cmake_minimum_required(VERSION 3.10)
project(MIS VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories(
    #VieCut/lib
    ${GUROBI_INCLUDE_DIR}
)

add_executable(MIS main.cpp algo/branch_and_reduce_algorithm.cpp algo/modified.cpp)
