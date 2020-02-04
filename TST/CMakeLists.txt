cmake_minimum_required(VERSION 3.15)
project(TST)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -v -lstdc++ -fopenmp -DNOCHANGE -DMKL_ILP64 -inline-force -qopt-report -xHost -O3 -Ofast -mkl")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -liomp5 -lpthread -ldl")


add_executable(TST main.cpp)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(TST Eigen3::Eigen)

find_package(CUDA REQUIRED)
target_include_directories(TST PUBLIC ${CUDA_INCLUDE_DIRS})
target_link_libraries(TST ${CUDA_LIBRARIES})
target_link_libraries(TST ${CUDA_cusolver_LIBRARY})
include_directories(${CUDA_INCLUDE_DIRS})
link_directories(${CUDA_LIBRARIES_DIRS})


find_package(MAGMA REQUIRED)
include_directories(${MAGMA_INCLUDE_DIRS})
link_directories(${MAGMA_LIBRARIES_DIRS})
target_include_directories(TST PUBLIC ${MAGMA_INCLUDE_DIRS})
target_link_libraries(TST ${MAGMA_LIBRARIES})

find_package(OpenMP REQUIRED)
target_link_libraries(TST OpenMP::OpenMP_CXX)


set(MKL_INCLUDE_DIR "/opt/intel/compilers_and_libraries_2019/linux/mkl/include")
set(MKL_LIBRARY_DIR "/opt/intel/compilers_and_libraries_2019/linux/mkl/lib")
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDE_DIR})
link_directories(${MKL_LIBRARY_DIR})
target_link_libraries(TST ${MKL_LIBRARY})


set(CMAKE_CXX_STANDARD 14)
