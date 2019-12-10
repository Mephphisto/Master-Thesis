//
//  main.cpp
//  TST
//
//  Created by Jakob Teuffel on 05.12.19.
//  Copyright © 2019 Jakob Teuffel. All rights reserved.
//

#define EIGEN_USE_MKL
#define OMP_NUM_THREADS 4
//#define DEBUG_ACTIVE
//#define EVAL_BY_CMD
#define EVAL_BY_CSV
constexpr auto MATRIX_SIZE = 200;
constexpr auto T_RES = 200;

#include <iostream>
#include <Eigen/Dense>
#include <omp.h>
#include <chrono>
#include <fstream>
#include <vector>
#include "_1DTopSuperConMatrix.hpp"



inline Eigen::Matrix3Xcd Build_Matrix()
{
    Eigen::MatrixXcd m(MATRIX_SIZE, MATRIX_SIZE);
#pragma omp parallel for
    for (auto i = 0; i < MATRIX_SIZE; i++) {
	//#pragma omp parallel for
        for (auto j = 0; j < MATRIX_SIZE; j++) {
            double j_d = double(j);
            double i_d = double(i);
            double z = j_d-i_d;
            double x = j_d + i_d;
            m(i, j) = std::complex<double>(1 + (x - MATRIX_SIZE / 2 - 2) / (z * z + 1), 0) * exp(std::complex<double>(0, z));
        }
    }
    return m;
}

int main()
{
#ifdef DEBUG_ACTIVE
	// get start time for runtime Calculation
	auto start = std::chrono::system_clock::now();

    // print number of used threads 
	std::cout << "Running " << Eigen::nbThreads() << " Eigen Threads" << std::endl;
#endif
	Eigen::MatrixXd All_EigenValues(MATRIX_SIZE, T_RES+1);
#pragma omp paralel num_threads(OMP_NUM_THREADS) private(k,m,Solver,EigenValues)
    {
        //Create Solver
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> Solver(MATRIX_SIZE);
        for (auto k = omp_get_thread_num(); k <= T_RES; k = OMP_NUM_THREADS){
            //Build Matrix to be Solved by MKL
            auto m = _1DTopSuperConMatrix(MATRIX_SIZE, -1 + 2* double(k) / T_RES, 0.1, 0.1).get();
            //Solve the Matrix
            Solver.compute(m);
            // Fetch  Eigenvalues from Solver and sort them in ascending order
            Eigen::VectorXd EigenValues = Solver.eigenvalues().col(0).real();
            std::sort(EigenValues.data(), EigenValues.data() + EigenValues.size());
            All_EigenValues.col(k) = EigenValues + Eigen::VectorXd::Constant(MATRIX_SIZE,(EigenValues.sum())/2);
        }
	}
#ifdef DEBUG_ACTIVE
	auto end = std::chrono::system_clock::now();
#endif

#ifdef EVAL_BY_CMD
	//print Eigenvalues to CMD
	std::cout << "eigenvalues:" << std::endl << EigenValues << std::endl;
#endif

#ifdef	EVAL_BY_CSV
	std::fstream csv_file("EigenValues_M"+ std::to_string(MATRIX_SIZE) + "_Tres"+ std::to_string(T_RES)+".csv",  std::fstream::out);
	assert(csv_file.is_open());
	csv_file << All_EigenValues;
	csv_file.close();
#endif
#ifdef DEBUG_ACTIVE
	std::cout << std::endl << std::endl << "Runtime = "
	<< std::chrono::duration_cast<std::chrono::milliseconds>(end- start).count()
	<< "ms" << std::endl;
#endif

#if defined EVAL_BY_CMD || defined DEBUG_ACTIVE
	//Wait for ENTER before Closing window
	std::cin.get();
#endif
}
