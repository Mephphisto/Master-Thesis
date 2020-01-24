cmake_minimum_required(VERSION 3.15)
project(TST)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -v -fopenmp -lstdc++")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ")

add_executable(TST main.cpp)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(TST Eigen3::Eigen)

set(CMAKE_CXX_STANDARD 14)
