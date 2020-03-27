cmake_minimum_required(VERSION 3.10)
project(TST)

set(CMAKE_THREAD_LIBS_INIT "-lpthread")
set(CMAKE_HAVE_THREADS_LIBRARY 1)
set(CMAKE_USE_WIN32_THREADS_INIT 0)
set(CMAKE_USE_PTHREADS_INIT 1)
set(THREADS_PREFER_PTHREAD_FLAG ON)

set(MKLROOT /opt/intel/mkl)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fopenmp -DMKL_ILP64 -m64 -I${MKLROOT}/include -DNOCHANGE")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl")
#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

set(CMAKE_CXX_STANDARD 14)

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

find_package(OpenMP)
if(OpenMP_FOUND)
    target_link_libraries(TST OpenMP::OpenMP_CXX)
endif()
