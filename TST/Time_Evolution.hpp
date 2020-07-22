//
// Created by Jakob Teuffel on 08.04.20.
//
#pragma once

#include <Eigen/Dense>
#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>
#include "Hamiltonian_Matrix.hpp"
#include "Parameters.hpp"
#include <cmath>
#include <omp.h>
#include "implicit_euler_Eigen.hpp"
#include <string>
/// Neccecary ? is this nneded for ODEint ?
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

// Eigen Typedefs
typedef Eigen::MatrixXcd Mat_cd;
typedef Eigen::VectorXcd Vec_cd;
typedef Eigen::VectorXd Vec;
typedef std::complex<double> cd;

//ODE Typedefs
typedef Vec_cd State;
typedef Mat_cd Jacobi;
typedef double Time;
typedef boost::numeric::odeint::implicit_euler_Eigen<cd> Stepper;

/// check is GPU is designated for diagonalisation
#if  defined(USE_MAGMA) && defined(USE_GPU)
#include "TST_Magma_Solver.hpp"
typedef magma::SelfAdjointEigenSolver<Mat_cd> MSolver;
#else
typedef Eigen::SelfAdjointEigenSolver<Mat_cd> MSolver;
#endif
// activate deactivate manual debugging diring execution
#undef DEBUG_ACTIVE2

/// Build QM Groundsate according to energy levels
inline Vec_cd Get_C_0(Vec eval, Mat_cd evec) {
    Vec_cd C_0 = Vec_cd(MATRIX_SIZE);
    // Summ up over Fermmie see and Majorana states thus including the first state above \mu
    for (size_t k = 0; k < MATRIX_SIZE; k++) {
        C_0 += (eval[k] < 10e-7) ? evec.col(k).normalized() : Vec_cd::Zero(MATRIX_SIZE);
    }
#ifdef DEBUG_ACTIVE2
    /// print out the Vectornorm thus the particle number in the Fermi see +1 (see above)
    std::cout << "|C_0| = " << C_0.norm() << std::endl;
#endif
    return C_0;
}

/// causes errors : Basis of compllex Vectorspace is overcomplete in Real subspace!
#deprecated
inline double Get_C_final_Overlap(Vec_cd const &C_f, Vec eval, Mat_cd evec) {
    double res = 0;
    // Summ up over Fermmie see and Majorana states
    for (size_t k = 0; k < MATRIX_SIZE; k++) {
        res += (eval[k] < 10e-7) ? std::abs(evec.col(k).dot(C_f)) : 0.0;
    }
    return res;
}

/// This computes the Energies from the Eigenvalues according to "Plegot & Ripka"
inline Vec Energys(Vec_cd eval, double tr) {
    Vec res = eval.col(0).real() - Vec::Constant(MATRIX_SIZE, 0.5 * (tr - eval.col(0).real().sum()));
    return res;
}


/// This is the ODE to be solved AKA the schroedinger Equation
template<typename T>
struct Schroedinger_of_cs {
private:
    double w;
public:
    explicit Schroedinger_of_cs(double omega) : w(omega) {}

    /// This is the ODE to be solved
    void operator()(const State &c, State &dcdt, Time theta) {
        dcdt = T(MATRIX_SIZE, T_COUPLE, theta , MU, DELTA).get() * c;
        dcdt *= cd(0, -1/w);
    }
};
/// This is the Jacobi Matrix of the Schroedinger Equation AKA the Hamiltonian. (Check this!)
template<typename T>
struct Schroedinger_of_cs_Jacobi {
private:
    double w;
public:
    explicit Schroedinger_of_cs_Jacobi(double omega) : w(omega) {}

    /// This is the ODE to be solved
    void operator()(const State &c, Jacobi &Jacobi, Time theta) {
        Jacobi = T(MATRIX_SIZE, T_COUPLE, theta, MU, DELTA).get();
    }
};

struct last_observer {
    State &m_state;
#ifdef DEBUG_ACTIVE2
    Vec_cd C_0;
#endif

    explicit last_observer(State &state) : m_state(state) {
#ifdef DEBUG_ACTIVE2
        C_0 = Vec_cd::Ones(MATRIX_SIZE);
#endif
    }
    void operator()(const State &x, double t) {
#ifdef DEBUG_ACTIVE2
        std::cout << "t= " << t << " |C_0-C_t| = " << (C_0 - x).norm() << " <C_0,C_t>/|C_0|² = "
                  << std::abs(C_0.dot(x)) / std::pow(C_0.norm(), 2) << std::endl;
        if (t - T_START < double(2.00) * M_PI / (T_RES + 1)) C_0 = x;
#endif
        std::string msg = "At k =  " + std::to_string(omp_get_thread_num()) +  "\t at t = " + std::to_string(t) + "\n";
        std::cout << msg;
        if (t == double(2) * M_PI + T_START) {
            m_state = x;
        }
    }
};

Vec Do_TE(Vec const &Omegas) {
    static_assert(std::is_base_of<Hamiltonian_Matrix, HAMILTONIAN>::value,
                  "Given Class is not derived from Hamiltonian_Matrix");
    Vec Rho_t(Omegas.size());
    Vec_cd C_0, eval;
    Mat_cd evec;
    {
        auto M = HAMILTONIAN(MATRIX_SIZE, T_COUPLE, T_START, MU, DELTA);
#ifdef DEBUG_ACTIVE2
        std::cout << "Hermiticity Check" << std::endl;
        M.verify_hermitiity();
        std::cout << "Matrix Det " << M.get().determinant() << std::endl;
        std::cout << "M*1= " << (M.get() * Vec_cd::Ones(MATRIX_SIZE)).norm() << std::endl;
#endif
        double tr = M.trace_A();
        MSolver Solver(MATRIX_SIZE);
        Solver.compute(M.get());
        assert(("Demoralisation :: No Success ", Solver.info() == Eigen::Success));
        eval = Solver.eigenvalues();
        evec = Solver.eigenvectors();
#ifdef DEBUG_ACTIVE2
        std::cout << " Debug Test Diagonalizeing Well?" << std::endl;
        Mat_cd D = (evec.adjoint() * M.get() * evec);
        D -= eval.asDiagonal();
        std::cout << "|THT^ -D| = " << D.norm() << std::endl;
        std::cout << "|(THT^ -D)*1| = " << (D * Vec_cd::Ones(MATRIX_SIZE)).norm() << std::endl;
        std::cout << "Debug end" << std::endl;
#endif

        C_0 = Get_C_0(Energys(eval, tr), evec);
    }
#pragma omp parallel for shared(Rho_t, C_0, eval, evec, Omegas, std::cout) default(none)
    for (size_t k = 0; k < Omegas.size(); k++) {
        Vec_cd C_End;
        Vec Enrgy_E;
        Mat_cd evec_E;
        {
            auto M = HAMILTONIAN(MATRIX_SIZE, T_COUPLE, T_END / Omegas[k], MU, DELTA);
            double tr_E = M.trace_A();
            MSolver Solver(MATRIX_SIZE);
            Solver.compute(M.get());
            if (Solver.info() != Eigen::Success) {
                std::cout << " Demoralisation :: No Success w= " << Omegas[k] << std::endl;
            }
            evec_E = Solver.eigenvectors();
            Enrgy_E = Energys(Solver.eigenvalues(), tr_E);
            C_End = Get_C_0(Enrgy_E, evec_E);
        }
        Vec_cd C_f(MATRIX_SIZE);
        State C = C_0;
        Schroedinger_of_cs<HAMILTONIAN> system(Omegas[k]);
        Schroedinger_of_cs_Jacobi<HAMILTONIAN> system2(Omegas[k]);
        Vec_cd dCdt(MATRIX_SIZE);
        boost::numeric::odeint::integrate_const(
                Stepper(),
                std::make_pair(system, system2),
                C,
                (double) T_START,
                (double) T_END,
                double(2.00) * M_PI / T_RES * Omegas[k],
                last_observer(C_f));
#ifdef DEBUG_ACTIVE2
        std::cout << "w= " << Omegas[k] << " |C_E-C_f| = " << (C_End - C_f).norm() << " <C_E,C_f>/|C_0|² = "
                  << std::abs(C_End.dot(C_f)) / std::pow(C_End.norm(), 2) << std::endl;
#endif
        Rho_t[k] = std::abs(C_0.dot(C_f)) / std::pow(C_0.norm(), 2);//Get_C_final_Overlap(C_f, Enrgy_E, evec_E);
    }
    return Rho_t;
}