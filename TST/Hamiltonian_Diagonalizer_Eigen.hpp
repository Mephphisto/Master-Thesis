//
//  Hamiltonian_Diagonalizer.hpp
//  TST
//
//  Created by Jakob Teuffel on 10.12.19.
//  Copyright © 2019 Jakob Teuffel. All rights reserved.
//

#pragma once

#include "Parameters.hpp"
#include "Hamiltonian_Matrix.hpp"
#include "Hamiltonian_Diagonalizer.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <chrono>
#include <fstream>
#include <vector>
#if OMP_NUM_THREADS > 1
#include <omp.h>
#endif


template<class T>
class Diagonalize_Hamiltonian_Eigen : public Diagonalize_Hamiltonian<T> {
     virtual void Compute() override {
#if OMP_NUM_THREADS > 1
#pragma omp parallel num_threads(OMP_NUM_THREADS)

        {
            // Get Thread ID
            int tid = omp_get_thread_num();
            // Create Solver
#else
            int tid = 0;
#endif
            Eigen::SelfAdjointEigenSolver<decltype(((T*)nullptr)->get())> Solver(MATRIX_SIZE);

            for (auto k = tid; k < T_RES; k += OMP_NUM_THREADS) {
#ifdef DEBUG_ACTIVE
                //verify accurate stride through sample space
                {
                    //Build String first and print then Prevents Race Condition!
                    std::string msg;
                    msg = "Thread " + std::to_string(tid) + " k= " + std::to_string(k) + "\n";
                    std::cout << msg;
                }
#endif
                //Storage for Matrix and Trace
                Eigen::MatrixXcd m;
                double tr;
                double t = T_START + (T_END - T_START) * (double(k) / double(T_RES));
                {
                    //Build Matrix to be Solved by MKL
                    T M = T(MATRIX_SIZE, T_COUPLE, t, MU, DELTA);
                    m = M.get();
                    tr = M.trace_A();
                }
                //Solve the Matrix
                Solver.compute(m);
                // Fetch  Eigenvalues from Solver
                this->All_EigenValues.col(k) = Solver.eigenvalues().col(0).real();
                this->t_s(k) = t;
                this->All_EigenVectors[k](Solver.eigenvectors());
            }
#if OMP_NUM_THREADS > 1
        }
#endif
    }
};