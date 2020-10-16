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
#include "Majorana_Transform.hpp"
#include "Typedefs.hpp"
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
                Mat_cd m;
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
                // Check successful Eigen solve
#pragma omp critical
                {
                    auto status = Solver.info();
                    if (status != Eigen::Success) {
                        if (status == Eigen::ComputationInfo::NoConvergence) std::cerr << "NoConvergence" << std::endl;
                        if (status == Eigen::ComputationInfo::InvalidInput) std::cerr << "InvalidInput"<< std::endl;
                        if (status == Eigen::ComputationInfo::NumericalIssue) std::cerr << "NumericalIssue"<< std::endl;
                    }
                }
                this->All_EigenValues.col(k) = Solver.eigenvalues().col(0).real();
                this->t_s(k) = t;
                this->All_EigenValues.col(k) = Solver.eigenvalues().col(0).real() +
                                               Vec::Constant(MATRIX_SIZE, (tr - Solver.eigenvalues().col(
                                                       0).real().sum()) / 2);

                size_t majoranas[2] = {0, 0};
                for (size_t t = 0; t < this->All_EigenValues.col(k).size(); t++) {
                    if (abs(this->All_EigenValues.col(k)[t]) < abs(this->All_EigenValues.col(k)[majoranas[0]])) {
                        majoranas[1] = majoranas[0];
                        majoranas[0] = t;
                    } else if (abs(this->All_EigenValues.col(k)[t]) < abs(this->All_EigenValues.col(k)[majoranas[1]])) {
                        majoranas[1] = t;
                    }
                }
#define MAJIZE
#ifdef MAJIZE
                auto tpl = Majoranaize(Solver.eigenvectors().col(majoranas[0]).normalized(),Solver.eigenvectors().col(majoranas[1]).normalized());
                this->All_EigenVectors.push_back(std::get<0>(tpl));
                this->All_EigenVectors.push_back(std::get<1>(tpl));
#else
                for (auto l : majoranas) {
                    std::cout << " EiVal[" << l << "] = " << this->All_EigenValues.col(k)[l];
                    this->All_EigenVectors.push_back(Solver.eigenvectors().col(l));
#endif
            }
#if OMP_NUM_THREADS > 1
        }
#endif
    }
};