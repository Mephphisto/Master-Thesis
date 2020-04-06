//
// Created by jakob on 02.03.20.
//

#pragma once

#include <math.h>

#define MLK_NUM_THREADS OMP_NUM_THREADS
//ICL NEEDS THIS
#define  _HAS_CONDITIONAL_EXPLICIT 0

//Choose a Solver

//#undef USE_OpenCL
//#define USE_GPU
#if MAGMA
#define USE_MAGMA
#endif
#ifdef USE_MAGMA
#define EIGEN_USE_LAPACKE_STRICT
#else
#define EIGEN_USE_MKL_ALL
#endif

//Define Data Aquisition
//#define DEBUG_ACTIVE
#define EVAL_EVAL_BY_CSV
#define EVAL_EVEC_BY_CSV
//#define SHOW_MATRIX_AND_QUIT
//#define MATRIX_TO_CSV
//#define MU_TO_CSV
//#define DELTA

//Deine Matrix Parameters

//#define PERIODIC_BOUNDRY
#define PHASE
//constexpr int GRID = 35;
//constexpr int  T_RES = 400;
constexpr int MATRIX_SIZE = (2 * GRID * GRID);
constexpr double T_START = 0.0000000001;
constexpr double T_END = 4.0 / static_cast<double>(GRID);
constexpr double T_COUPLE = -0.9;
constexpr double MU =  3;
constexpr double DELTA = 3;

