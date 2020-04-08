//
// Created by Jakob Teuffel on 08.04.20.
//
#pragma once

#include <Eigen/Dense>
#include "Hamiltonian_Matrix.hpp"
#include "Hamiltonian_MatrixSparse.hpp"
#include <omp.h>


typedef Eigen::MatrixXcd Mat_cd;
typedef Eigen::VectorXcd Vec_cd;
typedef Eigen::VectorXd Vec;
typedef std::complex<double> cd;

inline Vec_cd Get_State_0(Eigen::SelfAdjointEigenSolver<Mat_cd> Solver) {
    Vec_cd State = Eigen::VectroXcd.zeros(MATRIX_SIZE);
    auto eval = Solver.eigenvalues();
    auto evec = Solver.eigenvectors();
    // Summ up over Fermmie see and Majorana states
    for (size_t k = 0; k < MATRIX_SIZE; k++) {
        if (eval[k] < 10e-7) { //Numerical acuracy of Eigen Solvers
            State += evec.col(k);
        }
    }
    return State;
}

template<typename T>
inline Eigen::SelfAdjointEigenSolver<Mat_cd> Diagonalize(T Hamil) {
    bool right_T = std::is_base_of<Hamiltonian_Matrix, T>::value;
    right_T |= std::is_base_of<Hamiltonian_MatrixSparse, T>::value;
    assert(right_T);
    T M = T(MATRIX_SIZE, T_COUPLE, T_START, MU, DELTA);
    cd tr = M.trace_A();
    Eigen::SelfAdjointEigenSolver<Mat_cd> Solver(MATRIX_SIZE);
    Solver.compute(M.get());
    Solver.eigenvalues() = Solver.eigenvalues().col(0).real() +
                           Vec::Constant(MATRIX_SIZE, (tr - Solver.eigenvalues().col(
                                   0).real().sum()) / 2);
    return Solver;
}

inline Vec_cd Timeeval_fullrot(Vec_cd state_t, double w){

}

Vec Do(Vec Omegas){
    Vec_cd State_T0 = Get_State_0(Diagonalize<HAMILTONIAN>());
    Vec Rho_t(Omegas.size());
#pragma omp parallel for
    for (size_t k = 0 ; k < Omegas.size(); k++){
        Rho_t[k] = State_T0.adjoint()*Timeeval_fullrot(State_T0, Omegas[k]);
    }
    return Rho_t;
}