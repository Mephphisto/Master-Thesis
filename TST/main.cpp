//
//  main.cpp
//  TST
//
//  Created by Jakob Teuffel on 05.12.19.
//  Copyright © 2019 Jakob Teuffel. All rights reserved.
//

#define OMP_NUM_THREADS 16
//ICL NEEDS THIS
#define  _HAS_CONDITIONAL_EXPLICIT 0
#define DEBUG_ACTIVE
//#define EVAL_BY_CSV
//#define EVAL_EVAL_BY_CSV
#define EVAL_EVEC_BY_CSV
//#define SHOW_MATRIX_AND_QUIT
//#undef USE_OpenCL
#define USE_MAGMA
//#define USE_GPU
#define EIGEN_USE_LAPACKE
//#define MATRIX_TO_CSV
//#define PERIODIC_BOUNDRY
#define PHASE

constexpr auto MATRIX_SIZE = 2 * 100 * 100;
constexpr auto T_RES = 1;
constexpr double T_START = 1;
constexpr double T_END = 1;
constexpr double MU = -1;

constexpr double DELTA = 2;


#include "_1DTopSuperConMatrix.hpp"
#include "_2DTopSuperConMatrix.hpp"

#ifdef USE_MAGMA
#include "Hamiltonian_Diagonalizer_MAGMA.hpp"
#elif defined USE_OpenCL
#include "Hamiltonian_Diagonalizer_OpenCl.hpp"
#else
#include "Hamiltonian_Diagonalizer.hpp"
#endif

int main() {
#ifdef USE_MAGMA
    Diagonalize_Hamiltonian_magma<_2DTopSuperConMatrix>();
#elif defined USE_OpenCL
    Diagonalize_Hamiltonian_OpenCL<_2DTopSuperConMatrix>();
#else
    Diagonalize_Hamiltonian<_2DTopSuperConMatrix>();
#endif
    return 0;
}
