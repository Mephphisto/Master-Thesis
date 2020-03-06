//
// Created by jakob on 02.03.20.
//

#ifndef TST_PARAMETERS_HPP
#define TST_PARAMETERS_HPP

#include <math.h>

#define OMP_NUM_THREADS 16
//ICL NEEDS THIS
#define  _HAS_CONDITIONAL_EXPLICIT 0

//Choose a Solver

//#undef USE_OpenCL
//#define USE_MAGMA
//#define USE_GPU
//#define EIGEN_USE_LAPACKE

//Define Data Aquisition

#define DEBUG_ACTIVE
#define EVAL_EVAL_BY_CSV
#define EVAL_EVEC_BY_CSV
//#define SHOW_MATRIX_AND_QUIT
//#define MATRIX_TO_CSV
//#define MU_TO_CSV
//#define DELTA

//Deine Matrix Parameters

//#define PERIODIC_BOUNDRY
#define PHASE
constexpr int GRID = 32;
constexpr int MATRIX_SIZE = (2 * GRID * GRID);
constexpr int  T_RES = 3 * GRID;
constexpr double T_START = 0.0;
constexpr double T_END = ((double) T_RES + 1) / T_RES * M_PI;
constexpr double T_COUPLE = .95;
constexpr double MU =  -3;
constexpr double DELTA = 2;

#endif //TST_PARAMETERS_HPP
