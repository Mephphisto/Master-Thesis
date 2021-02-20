//
// Created by jakob on 02.03.20.
//

#pragma once

#include <cmath>

//Choose a Solver

//#undef USE_OpenCL
#if defined(USE_MAGMA) && !defined(USE_GPU)
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
//#define MAJIZE

//Deine Matrix Parameters
constexpr double pref = 1.0 / 3.0;
//#define PERIODIC_BOUNDRY
#define PHASE
constexpr int MATRIX_SIZE = (2 * GRID * GRID);
constexpr double T_START = -0.0;
constexpr double T_END = -(M_PI * 2.0 * T_ROT + T_START);
constexpr double T_COUPLE = pref * T_C;
constexpr double MU = pref * 3;
constexpr double DELTA = pref * 3;
constexpr double R_STIFF = 2;