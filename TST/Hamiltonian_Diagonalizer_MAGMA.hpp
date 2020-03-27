//
// Created by Jakob Teuffel on 18.12.19.
//

#pragma once
#define EIGEN_USE_LAPACKE

#include "Parameters.hpp"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#ifdef USE_GPU
#include "TST_Magma_Solver.hpp"
#endif
#include "Hamiltonian_Diagonalizer.hpp"
#include <magma_lapack.h>
#include <vector>
#include <chrono>


template<class T>
class Diagonalize_Hamiltonian_magma : public Diagonalize_Hamiltonian<T> {
    virtual void Compute() override {

#ifdef USE_GPU
        magma::SelfAdjointEigenSolver<Eigen::MatrixXcd> Solver(MATRIX_SIZE);
#else
        Eigen::SelfAdjointEigenSolver<decltype(((T *) nullptr)->get())> Solver(MATRIX_SIZE);
#endif
        for (auto k = 0; k < T_RES; k++) {
#ifdef DEBUG_ACTIVE
            //verify accurate stride through sample space
            {
                //Build String first and print then Prevents Race Condition!
                std::string msg;
                msg = " Status= " + std::to_string(static_cast<float>(k) / T_RES * 100) + "%\n";
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

            this->All_EigenValues.col(k) = Solver.eigenvalues().col(0).real() +
                                           Eigen::VectorXd::Constant(MATRIX_SIZE, (tr - Solver.eigenvalues().col(
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
            for (auto l : majoranas) {
                std::cout << " EiVal[" << l << "] = " << this->All_EigenValues.col(k)[l];
                this->All_EigenVectors.push_back(Solver.eigenvectors().col(l));
            }
            std::cout << std::endl;
            this->t_s(k) = t;
        }
    }

};