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
#include <boost/numeric/odeint/external/eigen/eigen.hpp>

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
typedef boost::numeric::odeint::runge_kutta4<Vec_cd> stepper_type;


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
inline std::tuple<Vec, Mat_cd> Diagonalize() {
    static_assert(std::is_base_of<Hamiltonian_Matrix, T>::value);
    T M = T(MATRIX_SIZE, T_COUPLE, T_START, MU, DELTA);
    double tr = M.trace_A();
    MSolver Solver(MATRIX_SIZE);
    Solver.compute(M.get());
    Vec eval = Solver.eigenvalues().col(0).real() +
               Vec::Constant(MATRIX_SIZE, 0.5 * (tr - (Solver.eigenvalues().col(0).real()).sum()));;
    Mat_cd evec = Solver.eigenvectors();
    return std::make_tuple(eval, evec);
}

/// This is the ODE to be solved
template<typename T>
struct Schroedinger_of_cs {
    static_assert(std::is_base_of<Hamiltonian_Matrix, T>::value);
private:
    Mat_cd evec;
    double w;
public:
    Schroedinger_of_cs( double omega, Mat_cd eigen_vectors): w(omega), evec(eigen_vectors) {}

    /// This is the ODE to be solved
    void operator()(Vec_cd &c, Vec_cd &dcdt, double theta) {
        auto get = T(MATRIX_SIZE, T_COUPLE, theta, MU, DELTA).get();

        Vec_cd vm;
        dcdt = Vec_cd::Zero(MATRIX_SIZE);
        for (size_t l = 0; l < MATRIX_SIZE; l++) {
            vm.noalias() = evec.col(l).adjoint() * get;
            for (size_t k = 0; k < MATRIX_SIZE; k++) {
                dcdt[l] += vm.dot(evec.col(k)) * c[k];
            }
        }
        dcdt *= cd(0.0, -1.0 / w);


         //dcdt = cd(0, -1.0 / w) * evec.conjugate() * get * evec * c;
    }
};


struct last_observer {
    Vec_cd  &m_state;

    last_observer( Vec_cd &state ) : m_state(state){
    }

    void operator()( const Vec_cd &x , double t )
    {
        m_state = x ;
    }
};

Vec Do(Vec Omegas) {
    Vec eval;
    Mat_cd evec;
    std::tie(eval, evec) = Diagonalize<HAMILTONIAN>();
    Vec_cd C_0 = Get_C_0(eval);
    Vec Rho_t(Omegas.size());
#pragma omp parallel for
    for (size_t k = 0; k < Omegas.size(); k++) {
        Vec_cd C_f(MATRIX_SIZE);
        Schroedinger_of_cs<HAMILTONIAN> system( Omegas[k], evec);
        boost::numeric::odeint::integrate_const(
                stepper_type(),
                system,
                C_0,
                double(0),
                double(2) * M_PI,
                double(2.00) * M_PI / T_RES,
                last_observer(C_f));

        cd aux = (C_0.adjoint() * C_f).value();
        Rho_t[k] = std::abs(aux);
    }
    return Rho_t;
}