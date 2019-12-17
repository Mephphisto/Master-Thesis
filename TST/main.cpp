//
//  main.cpp
//  TST
//
//  Created by Jakob Teuffel on 05.12.19.
//  Copyright © 2019 Jakob Teuffel. All rights reserved.
//


#define EIGEN_USE_MKL
#define OMP_NUM_THREADS 16
//ICL NEEDS THIS
#define  _HAS_CONDITIONAL_EXPLICIT 0
//#define DEBUG_ACTIVE
//#define EVAL_BY_CMD
#define EVAL_BY_CSV
//#define SHOW_MATRIX_AND_QUIT

constexpr auto MATRIX_SIZE = 400;
constexpr auto T_RES = 64;
constexpr double T_START = 0;
constexpr double T_END = 6;
constexpr double MU = 0;
constexpr double DELTA = 1;

#include "_1DTopSuperConMatrix.hpp"
#include "Hamiltonian_Diagonalizer.hpp"



int main()
{
    Diagonalize_Hamiltonian<_1DTopSuperConMatrix>();
    return 0;
}
