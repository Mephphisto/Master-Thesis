//
// Created by Jakob Teuffel on 08.04.20.
//
#pragma once

#include <Eigen/Dense>
#include <omp.h>
#include "Hamiltonian_Matrix.hpp"
#include "Hamiltonian_MatrixSparse.hpp"
#include "Parameters.hpp"
#include "_1DTopSuperConMatrix.hpp"
#include "_2DTopSuperConMatrix.hpp"
#include "_2DTopSuperConMatrixSparse.hpp"
#include "Time_Evolution.hpp"
#include <boost/numeric/odeint.hpp>
#include <math.h>

#if  defined(USE_MAGMA) && defined(USE_GPU)

#include "TST_Magma_Solver.hpp"

#endif
typedef Eigen::MatrixXcd Mat_cd;
typedef Eigen::VectorXcd Vec_cd;
typedef Eigen::VectorXd Vec;
typedef std::complex<double> cd;
#if  defined(USE_MAGMA) && defined(USE_GPU)
typedef magma::SelfAdjointEigenSolver<Mat_cd> MSolver;
#else
typedef Eigen::SelfAdjointEigenSolver<Mat_cd> MSolver;
#endif


namespace pl = std::placeholders;

inline Vec_cd Get_C_0(Vec eval) {
    Vec_cd C_0 = Vec_cd(MATRIX_SIZE);
    // Summ up over Fermmie see and Majorana states
    for (size_t k = 0; k < MATRIX_SIZE; k++) {
        C_0[k] = (eval[k] < 10e-7) ? cd(1.0) : cd(0.0);
    }
    return C_0;
}

template<typename T>
inline auto Diagonalize() {
    static_assert(std::is_base_of<Hamiltonian_Matrix, T>::value);
    T M = T(MATRIX_SIZE, T_COUPLE, T_START, MU, DELTA);
    double tr = M.trace_A();
    MSolver Solver(MATRIX_SIZE);
    Solver.compute(M.get());
    return std::tie(Solver, tr);
}

/// This is the ODE to be solved
template<typename T>
class Solve_Timeeval {
    static_assert(std::is_base_of<Hamiltonian_Matrix, T>::value);
private:
    Vec_cd evec, initial_c;
    double w;
public:
    Solve_Timeeval(Vec_cd initial, double omega, Vec_cd eigen_vectors) {
        this->initial_c = initial;
        this->evec = eigen_vectors;
        this->w = omega;
    }

    /// This is the ODE to be solved
    void Schroedinger(Vec_cd c, Vec_cd dcdt, double theta) {

        dcdt = Vec_cd::Zero(MATRIX_SIZE);
        for (size_t l = 0; l < MATRIX_SIZE; l++) {
            for (size_t k = 0; k < MATRIX_SIZE; k++) {
                dcdt[k] += (evec.col(k).adjoint() * T(MATRIX_SIZE, T_COUPLE, theta, MU, DELTA).get() * evec.col(l)).value() * c[l];
            }
        }
        dcdt /= w;
    }
};


Vec Do(Vec Omegas) {
    auto PSolver = Diagonalize<HAMILTONIAN>();
    auto Solver = std::get<0>(PSolver);
    double tr = std::get<1>(PSolver);
    Vec eval = Solver.eigenvalues().col(0).real() +
               Vec::Constant(MATRIX_SIZE, 0.5 * (tr - (Solver.eigenvalues().col(0).real()).sum()));;
    auto evec = Solver.eigenvectors();
    Vec_cd C_0 = Get_C_0(eval);
    Vec Rho_t(Omegas.size());
#pragma omp parallel for
    for (size_t k = 0; k < Omegas.size(); k++) {
        Vec_cd state_c = C_0;
        boost::numeric::odeint::integrate(std::bind( &Solve_Timeeval<HAMILTONIAN>::Schroedinger , Solve_Timeeval<HAMILTONIAN>(C_0, Omegas[k], evec) , pl::_1 , pl::_2 , pl::_3 )
                , state_c, double(0), double(2) * M_PI, double(2.00) * M_PI /T_RES);
        cd aux = (C_0.adjoint() * state_c).value();
        Rho_t[k] = std::abs(aux);
    }
    return Rho_t;
}