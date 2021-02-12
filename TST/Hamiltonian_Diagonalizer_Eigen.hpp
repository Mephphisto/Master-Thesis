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
#ifdef MAJIZE
#include "Majoranize.hpp"
#endif

#if OMP_NUM_THREADS > 1

#include <omp.h>
#include "Majorana_Transform.hpp"
#include "Typedefs.hpp"

#endif


template<class T>
class Diagonalize_Hamiltonian_Eigen : public Diagonalize_Hamiltonian<T> {
    inline void Compute() override {
#if OMP_NUM_THREADS > 1
#pragma omp parallel num_threads(OMP_NUM_THREADS)

        {
            // Get Thread ID
            int tid = omp_get_thread_num();
            // Create Solver
#else
            int tid = 0;
#endif
            Eigen::SelfAdjointEigenSolver<Mat_cd> Solver(MATRIX_SIZE);


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
                Mat_cd m;
                double tr;
                double t = T_START + (T_END - T_START) * ((k == 0) ? 0.0 : (double(k) / double(T_RES - 1)));
                {
                    //Build Matrix to be Solved by MKL
                    T M = T(MATRIX_SIZE, T_COUPLE, t, MU, DELTA);
                    m = M.get();
                    tr = 0.0;//M.trace_A();
                }
                //Solve the Matrix
                Solver.compute(m);
                // Fetch  Eigenvalues from Solver
                // Check successful Eigen solve
#pragma omp critical
                {
                    if (Solver.info() != Eigen::Success) {
                        if (Solver.info() == Eigen::ComputationInfo::NoConvergence)
                            std::cerr << "NoConvergence" << std::endl;
                        if (Solver.info() == Eigen::ComputationInfo::InvalidInput)
                            std::cerr << "InvalidInput" << std::endl;
                        if (Solver.info() == Eigen::ComputationInfo::NumericalIssue)
                            std::cerr << "NumericalIssue" << std::endl;
                    }
                }
                this->All_EigenValues.col(k) = Solver.eigenvalues().col(0).real();
                this->t_s(k) = t;
                this->All_EigenValues.col(k) = Solver.eigenvalues().col(0).real() +
                                               Vec::Constant(MATRIX_SIZE, (tr - Solver.eigenvalues().col(
                                                       0).real().sum()) / 2);

                size_t majoranas[2] = {0, 0};

#pragma unroll
                for (size_t l = 0; l < this->All_EigenValues.col(k).size(); l++) {
                    if (abs(this->All_EigenValues.col(k)[l]) < abs(this->All_EigenValues.col(k)[majoranas[0]])) {
                        majoranas[1] = majoranas[0];
                        majoranas[0] = l;
                    } else if (abs(this->All_EigenValues.col(k)[l]) <
                               abs(this->All_EigenValues.col(k)[majoranas[1]])) {
                        majoranas[1] = l;
                    }
                }
#ifdef MAJIZE
                auto tpl = Fermiize(Solver.eigenvectors().col(majoranas[0]).normalized(),
                                    Solver.eigenvectors().col(majoranas[1]).normalized());
                this->All_EigenVectors.push_back(std::get<0>(tpl));
                this->All_EigenVectors.push_back(std::get<1>(tpl));
#else

#pragma unroll
                for (auto l : majoranas) {
#pragma omp critical
                    std::cout << " EiVal[" << l << "] = " << this->All_EigenValues.col(k)[l];
                    this->All_EigenVectors.push_back(Solver.eigenvectors().col(l));
                    }
#endif
            }
#if OMP_NUM_THREADS > 1
        }
#endif
    }
};