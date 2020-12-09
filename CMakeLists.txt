cmake_minimum_required(VERSION 3.10)
project(MIS VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories(
    Metis/build/Linux-x86_64/libmetis
)

#add_subdirectory(InertialFlowCutter)

add_executable(MIS main.cpp algo/branch_and_reduce_algorithm.cpp algo/modified.cpp)

target_link_libraries(MIS ${CMAKE_SOURCE_DIR}/Metis/build/Linux-x86_64/libmetis/libmetis.a)
#target_link_libraries(MIS inertial_flow)

