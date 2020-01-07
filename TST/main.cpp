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
//#define DEBUG_ACTIVE
//#define EVAL_BY_CMD
#define EVAL_BY_CSV
#define SHOW_MATRIX_AND_QUIT

constexpr auto MATRIX_SIZE = 12;
constexpr auto T_RES = 64;
constexpr double T_START = 3;
constexpr double T_END = 6;
constexpr double MU = 2;
constexpr double DELTA = 1;

#include "_1DTopSuperConMatrix.hpp"
#include "_2DTopSuperConMatrix.hpp"
#include "Hamiltonian_Diagonalizer.hpp"
#include "Hamiltonian_Diagonalizer_OpenCl.hpp"



int main()
{
#ifdef USE_OpenCL
    Diagonalize_Hamiltonian_OpenCL<_2DTopSuperConMatrix>();
#else
    //Diagonalize_Hamiltonian<_2DTopSuperConMatrix>();
#endif
    return 0;
}
