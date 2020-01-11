cmake_minimum_required(VERSION 3.15)
project(TST)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -v -lstdc++ -fopenmp")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -v ")
include_directories("/opt/clBLAS/include/" "/opt/clBLAS/bin/" "/opt/clBLAS/lib/")
LINK_DIRECTORIES("/opt/clBLAS/include/" "/opt/clBLAS/bin/" "/opt/clBLAS/lib/")

set(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} /opt/clBLAS/lib/)

add_executable(TST main.cpp)

option(EUSE_OCL "Decide weather or not to use OpenCL" ON)
IF(EUSE_OCL)
    ADD_DEFINITIONS(-DUSE_OpenCL)
    find_library( clBLAS REQUIRED)
    target_link_libraries(TST clBLAS)
    find_package( OpenCL ${OPENCL_VERSION} )
    target_link_libraries(TST OpenCL::OpenCL)
ENDIF(EUSE_OCL)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries (TST Eigen3::Eigen)

option(EUSE_MKL "Decide weather or not to use MKL" OFF)
IF(EUSE_MKL)
    find_package(BLAS REQUIRED)
    ADD_DEFINITIONS(-DEIGEN_USE_MKL)
ENDIF(EUSE_MKL)



set(CMAKE_CXX_STANDARD 14)