cmake_minimum_required(VERSION 3.10)
project(TST)

if (FALSE)
    set(CMAKE_THREAD_LIBS_INIT "-pthread")
    set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_USE_WIN32_THREADS_INIT 0)
    set(CMAKE_USE_PTHREADS_INIT 1)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
endif ()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(ENV{MKLROOT} "/opt/intel/compilers_and_libraries_2019.5.281/linux/mkl")
set(ENV{LD_LIBRARY_PATH} "/opt/intel/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64_lin:/opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/libfabric/lib:/opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/lib/release:/opt/intel/compilers_and_libraries_2019.5.281/linux/mpi/intel64/lib:/opt/intel/compilers_and_libraries_2019.5.281/linux/ipp/lib/intel64:/opt/intel/compilers_and_libraries_2019.5.281/linux/compiler/lib/intel64_lin:/opt/intel/compilers_and_libraries_2019.5.281/linux/mkl/lib/intel64_lin:/opt/intel/compilers_and_libraries_2019.5.281/linux/tbb/lib/intel64/gcc4.7:/opt/intel/compilers_and_libraries_2019.5.281/linux/tbb/lib/intel64/gcc4.7:/opt/intel/debugger_2019/libipt/intel64/lib:/opt/intel/compilers_and_libraries_2019.5.281/linux/daal/lib/intel64_lin:/opt/intel/compilers_and_libraries_2019.5.281/linux/daal/../tbb/lib/intel64_lin/gcc4.4")
set(ENV{MKL_NUM_THREADS} ${THREADS})
set(ENV{OMP_NUM_THREADS 1) #${THREADS})
set(ENV{VECLIB_MAXIMUM_THREADS} ${THREADS})



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -DNOCHANGE")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ldl")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")

set(CMAKE_CXX_STANDARD 14)

add_executable(TST main.cpp)

target_compile_definitions(TST PUBLIC OMP_NUM_THREADS=${THREADS})
target_compile_definitions(TST PUBLIC GRID=${GRID})
target_compile_definitions(TST PUBLIC T_RES=${T_RES})
target_compile_definitions(TST PUBLIC MKL_NUM_THREADS=${THREADS})
if (USE_MAGMA)
    target_compile_definitions(TST PUBLIC USE_MAGMA="")
endif ()
if (USE_GPU)
    target_compile_definitions(TST PUBLIC USE_GPU="")
endif ()
if(VERBOSE)
target_compile_definitions(TST PUBLIC DEBUG_ACTIVE="")
endif()
if(TIME_EVOLUTION)
    target_compile_definitions(TST PUBLIC TIME_EVOLUTION="")
endif()


find_package(Eigen3 3.3 REQUIRED NO_MODULE)
target_link_libraries(TST Eigen3::Eigen)
if (${USE_MAGMA})
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
endif ()

find_package(MKL REQUIRED)
if (MKL_FOUND)
    include_directories(${MKL_INCLUDE_DIR})
    target_link_libraries(TST ${MKL_LIBRARIES} ${MKL_OPENMP_LIBRARY})
endif ()

find_package(OpenMP)
if (OPENMP_FOUND)
    target_link_libraries(TST OpenMP::OpenMP_CXX)
    link_directories(${OpenMP_iomp5_LIBRARY} ${OpenMP_CXX_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()


find_package(Boost 1.72.0 REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
    message(${Boost_INCLUDE_DIRS})
endif()
