//
//  Defs.hpp
//  sim
//
//  Created by Jakob Teuffel on 10.12.19.
//  Copyright © 2019 Jakob Teuffel. All rights reserved.
//

#pragma once
#include <cmath>
//Choose a Solver
// USE MKL Beneath eigen
#define EIGEN_USE_MKL_ALL
// Choose HAMILTONIAN MATRIX CLASS
#define HAMILTONIAN _2DTopSuperConMatrix

//Define Data Aquisition for diagonalize by defining undefining
// Save Eigen Values
//#define EVAL_EVAL_BY_CSV

// Save lowest Eigen vectors
//#define EVAL_EVEC_BY_CSV

// Save Matrix entries
//#define MATRIX_TO_CSV

// Save chemical Potential
//#define MU_TO_CSV

// Save chemical Pair creation term
//#define DELTA_TO_CSV


// Debug Options
//#define SHOW_MATRIX_AND_QUIT


//Deine Matrix Parameters
//Vortex Radius
constexpr double RADIUS = GRID * 1.0/8.0;

//Use Periodic Boundary's
//#define PERIODIC_BOUNDRY

// Use Phasse winding
#define PHASE

// use Gaussian Vortices
#define GAUSS

// Define Shift
#define SHIFT 4.0 *t
//Grid size
constexpr int GRID = 35;
//Matrix size
constexpr int MATRIX_SIZE = (2 * GRID * GRID);
// start angle
constexpr double T_START = -0.0;
// % of full roatation
constexpr double ROT = 0.5;
// End time
constexpr double T_END = -(M_PI * 2.0 * T_ROT + T_START);
// hopping term
constexpr double T_COUPLE = 0.5;
// Chemical potential Mu_0
constexpr double MU =  3.0;
// Pair creation term
constexpr double DELTA =   1.0;
// number of time steps / THREADSHL
constexpr int T_RES=4

// type defs to simplyfy matrix definition
typedef Eigen::MatrixXcd Mat_cd;
typedef Eigen::MatrixXd Mat;
typedef Eigen::VectorXcd Vec_cd;
typedef Eigen::VectorXd Vec;
typedef std::complex<double> cd;

// Low level Mkl threads 
#define THREADS 4
// High Level OpenMP threads
#define THREADSHL 2
