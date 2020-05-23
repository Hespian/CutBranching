cmake_minimum_required(VERSION 3.10)
project(MIS VERSION 1.0.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include_directories(
    #VieCut/lib
    Metis/build/Linux-x86_64/libmetis
)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/NestedDissection/lib)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/NestedDissection/lib/io)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/NestedDissection/lib/partition)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/NestedDissection/lib/partition/uncoarsening/refinement/quotient_graph_refinement/flow_refinement)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/NestedDissection/lib/tools)

add_executable(MIS main.cpp algo/branch_and_reduce_algorithm.cpp algo/modified.cpp NestedDissection/lib/dissection/reductions.cpp 
NestedDissection/lib/dissection/ordering_tools.cpp NestedDissection/extern/argtable3-3.0.3/argtable3.c)
#add_executable(MIS oct_bench.cpp algo/branch_and_reduce_algorithm.cpp algo/modified.cpp)
#add_executable(MIS pg2.cpp algo/branch_and_reduce_algorithm.cpp algo/modified.cpp RND/lib/dissection/reductions.cpp 
#               RND/lib/dissection/ordering_tools.cpp RND/extern/argtable3-3.0.3/argtable3.c)
target_link_libraries(MIS ${CMAKE_SOURCE_DIR}/Metis/build/Linux-x86_64/libmetis/libmetis.a)

