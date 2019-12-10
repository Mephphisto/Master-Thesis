//
//  main.cpp
//  TST
//
//  Created by Jakob Teuffel on 05.12.19.
//  Copyright © 2019 Jakob Teuffel. All rights reserved.
//

#define EIGEN_USE_MKL
#define OMP_NUM_THREADS 16
//#define DEBUG_ACTIVE
//#define EVAL_BY_CMD
#define EVAL_BY_CSV
#define SHOW_MATRIX_AND_QUIT

constexpr auto MATRIX_SIZE = 10;
constexpr auto T_RES = 2 * OMP_NUM_THREADS;
constexpr double T_START = -0.6;
constexpr double T_END = 0.6;
constexpr double MU = 0.1;
constexpr double DELTA = 0.1;

#include <iostream>
#include <Eigen/Dense>
#include <omp.h>
#include <chrono>
#include <fstream>
#include "_1DTopSuperConMatrix.hpp"
#include "Hamiltonian_Diagonalizer.hpp"



int main()
{
    Diagonalize_Hamiltonian<_1DTopSuperConMatrix>();
    return 0;
}
