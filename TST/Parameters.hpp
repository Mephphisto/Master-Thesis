//
// Created by jakob on 02.03.20.
//

#pragma once

#include <math.h>

#define OMP_NUM_THREADS 16
#define MLK_NUM_THREADS OMP_NUM_THREADS
//ICL NEEDS THIS
#define  _HAS_CONDITIONAL_EXPLICIT 0

//Choose a Solver

//#undef USE_OpenCL
#define USE_MAGMA
//#define USE_GPU

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
//constexpr int GRID = 32;
constexpr int MATRIX_SIZE = 2*101*100;// (2 * GRID * GRID);
constexpr int  T_RES = 400; //3 * GRID;
constexpr double T_START = -1.0;
constexpr double T_END = 1.0;
constexpr double T_COUPLE = 3.0;
constexpr double MU =  1;
constexpr double DELTA = 2;

