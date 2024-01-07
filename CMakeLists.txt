cmake_minimum_required(VERSION 3.12.0)
project(NMPDE-Project LANGUAGES CXX C)

include(common/cmake-common.cmake)

add_executable(Project_Parallel src/Project_Parallel.cpp src/Poisson3D_Parallel.cpp)
deal_ii_setup_target(Project_Parallel)
