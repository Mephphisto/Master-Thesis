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
#define EVAL_BY_CSV
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
//#define MU_TO_CSV
//#define DELTA_TO_CSV

#include <math.h>

constexpr auto MATRIX_SIZE = 2 * 85 * 85;
constexpr auto T_RES = 32;
constexpr double T_START = 0.01;
constexpr double T_END = ((double)T_RES+1)/T_RES * M_PI;
constexpr double T_COUPLE = 1;
constexpr double MU = 2;
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
