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


/// This is the ODE to be solved
template<typename T>
struct Schroedinger_of_cs {
    static_assert(std::is_base_of<Hamiltonian_Matrix, T>::value);
private:
    Mat_cd evec;
    double w, Norm;
public:
    Schroedinger_of_cs(double omega, Mat_cd eigen_vectors, double n) : w(omega), evec(eigen_vectors), Norm(n) {}

    /// This is the ODE to be solved
    void operator()(Vec_cd &c, Vec_cd &dcdt, double theta) {
        auto H = T(MATRIX_SIZE, T_COUPLE, theta, MU, DELTA).get()/Norm;
        dcdt = cd(0, -1 / w) * evec.adjoint() * H * evec * c;
    }
};

struct last_observer {
    Vec_cd &m_state;

    last_observer(Vec_cd &state) : m_state(state) {
    }

    void operator()(const Vec_cd &x, double t) {
        m_state = x;
    }
};

Vec Do(Vec Omegas) {
    static_assert(std::is_base_of<Hamiltonian_Matrix, HAMILTONIAN>::value);
    Vec Rho_t(Omegas.size());
    Vec_cd C_0, eval;
    double Norm;
    Mat_cd evec;
    {
        auto M = HAMILTONIAN(MATRIX_SIZE, T_COUPLE, T_START, MU, DELTA);
        double tr = M.trace_A();
        cd det = M.get().determinant();
        Norm = std::abs(pow(det,1/MATRIX_SIZE));
        //MSolver Solver(MATRIX_SIZE);
        //Solver.compute(M.get()/Norm);
#ifdef DEBUG_ACTIVE
        std::cout << "det = " << det << std::endl
        << "Norm = " << Norm << std::endl
        <<"det(m) = " << (M.get()/Norm).determinant() << std::endl;
#endif
/*        std::cout << " Solver Success? " << (Solver.info() == Eigen::NoConvergence) << std::endl;
        eval = Solver.eigenvalues();
        evec = Solver.eigenvectors();
        C_0 = Get_C_0(eval.col(0).real() -
                      Vec::Constant(MATRIX_SIZE, 0.5 * (tr - (Solver.eigenvalues().col(0).real()).sum())));
    */}/*

#pragma omp parallel for
    for (size_t k = 0; k < Omegas.size(); k++) {
        Vec_cd C_f(MATRIX_SIZE);
        Schroedinger_of_cs<HAMILTONIAN> system(Omegas[k], evec, Norm);
        boost::numeric::odeint::integrate_const(
                stepper_type(),
                system,
                C_0,
                T_START,
                double(2) * M_PI + T_START,
                double(2.00) * M_PI / T_RES,
                last_observer(C_f));
        Rho_t[k] = std::abs(C_0.normalized().dot(C_f.normalized())) ;
    }*/
    return Rho_t;
}