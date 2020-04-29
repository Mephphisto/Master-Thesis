//
// Created by jakob on 02.03.20.
//

#pragma once

#include <math.h>

//ICL NEEDS THIS
#define  _HAS_CONDITIONAL_EXPLICIT 0

//Choose a Solver

//#undef USE_OpenCL
#if defined(USE_MAGMA)
#define EIGEN_USE_LAPACKE_STRICT
#else
#define EIGEN_USE_MKL_ALL
#endif

#define HAMILTONIAN _2DTopSuperConMatrix

//Define Data Aquisition
#define EVAL_EVAL_BY_CSV
#define EVAL_EVEC_BY_CSV
//#define SHOW_MATRIX_AND_QUIT
//#define MATRIX_TO_CSV
//#define MU_TO_CSV
//#define DELTA_To_CSV

//Deine Matrix Parameters
constexpr double pref  = 0.4;
//#define PERIODIC_BOUNDRY
#define PHASE
//constexpr int GRID = 35;
//constexpr int  T_RES = 400;
constexpr int MATRIX_SIZE = (2 * GRID * GRID);
constexpr double T_START = 0.2;
constexpr double T_END = 1.0 / static_cast<double>(GRID);
constexpr double T_COUPLE = -0.9*pref;
constexpr double MU =  pref*3;
constexpr double DELTA = pref*3;

