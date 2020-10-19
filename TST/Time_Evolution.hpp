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
#include "Majorana_Transform.hpp"
#include "Typedefs.hpp"
#include <cmath>
#include <omp.h>
//#include "implicit_euler_Eigen.hpp"
#include <string>

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>


typedef double Time;


#if  defined(USE_MAGMA) && defined(USE_GPU)
#include "TST_Magma_Solver.hpp"
typedef magma::SelfAdjointEigenSolver<Mat_cd> MSolver;
#else
typedef Eigen::SelfAdjointEigenSolver<Mat_cd> MSolver;
#endif

inline std::tuple<Vec_cd, Vec_cd, Vec_cd> Get_C_0_Majs(const Vec &eval, const Mat_cd &evec) {
    Vec_cd C_0 = Vec_cd(MATRIX_SIZE);
    // Summ up over Fermmie see and Majorana states
    size_t majoranas[2] = {0, 0};

#pragma unroll
    for (size_t t = 0; t < eval.size(); t++) {
        if (abs(eval[t]) < abs(eval[majoranas[0]])) {
            majoranas[1] = majoranas[0];
            majoranas[0] = t;
        } else if (abs(eval[t]) < abs(eval[majoranas[1]])) {
            majoranas[1] = t;
        }
    }

#pragma unroll
    for (size_t k = 0; k < MATRIX_SIZE; k++) {
        if ((k != majoranas[0]) || k > (k != majoranas[1])) {
            C_0 += (eval[k] < 10e-7) ? evec.col(k).normalized() : Vec_cd::Zero(MATRIX_SIZE);
        }
    }
#ifdef DEBUG_ACTIVE
    std::cout << "E_maj1 = " << eval[majoranas[0]] << " E_maj2 = " << eval[majoranas[1]] << std::endl;
#endif
    auto tpl = Majoranaize(evec.col(majoranas[0]).normalized(), evec.col(majoranas[1]).normalized());
    C_0 += std::get<0>(tpl) + std::get<1>(tpl);
#ifdef DEBUG_ACTIVE
    std::cout << "|C_0| = " << C_0.norm() << std::endl;
#endif
    return std::make_tuple(C_0, std::get<0>(tpl), std::get<1>(tpl));
}

inline double Get_C_final_Overlap(Vec_cd const &C_f, const Vec &eval, const Mat_cd &evec) {
    double res = 0;
    // Summ up over Fermmie see and Majorana states

#pragma unroll
    for (size_t k = 0; k < MATRIX_SIZE; k++) {
        res += (eval[k] < 10e-7) ? std::abs(evec.col(k).dot(C_f)) : 0.0;
    }
    return res;
}

inline Vec Energys(const Vec_cd &eval, const double &tr) {
    Vec res = eval.col(0).real() - Vec::Constant(MATRIX_SIZE, 0.5 * (tr - eval.col(0).real().sum()));
    return res;
}


/// This is the ODE to be solved
template<typename T>
struct Schroedinger_of_cs {
private:
    double w;
public:
    explicit Schroedinger_of_cs(const double &omega) : w(omega) {}

    /// This is the ODE to be solved
    inline void operator()(const Vec_cd &c, Vec_cd &dcdt, const Time theta) {
        dcdt = T(MATRIX_SIZE, T_COUPLE, theta, MU, DELTA).get() * c;
        dcdt *= cd(0, -1 / w);
    }
};

/*
template<typename T>
struct Schroedinger_of_cs_Jacobi {
private:
    double w;
public:
    explicit Schroedinger_of_cs_Jacobi(double omega) : w(omega) {}

    /// This is the ODE to be solved
    inline void operator()(const State &c, Jacobi &Jacobi, Time theta) {Vec_cd
        auto M (T(MATRIX_SIZE, T_COUPLE, theta, MU, DELTA).get());
        Jacobi = T(MATRIX_SIZE, T_COUPLE, theta, MU, DELTA).get();
    }
};
*/
struct last_observer {
    Vec_cd &m_state;

    explicit last_observer(Vec_cd &state) : m_state(state) {
    }

    void operator()(const Vec_cd &x, const double &t) {
        m_state = x;
    }
};

Mat Do_TE(Vec const &Omegas) {
    static_assert(std::is_base_of<Hamiltonian_Matrix, HAMILTONIAN>::value,
                  "Given Class is not derived from Hamiltonian_Matrix");
    Mat Rho_t(3, Omegas.size());
    Vec_cd C_0, eval, Maj1, Maj2;
    Mat_cd evec;
    {
        auto M = HAMILTONIAN(MATRIX_SIZE, T_COUPLE, T_START, MU, DELTA);
#ifdef DEBUG_ACTIVE
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
#ifdef DEBUG_ACTIVE
        std::cout << " Debug Test Diagonalizeing Well?" << std::endl;
        Mat_cd D = (evec.adjoint() * M.get() * evec);
        D -= eval.asDiagonal();
        std::cout << "|THT^ -D| = " << D.norm() << std::endl;
        std::cout << "|(THT^ -D)*1| = " << (D * Vec_cd::Ones(MATRIX_SIZE)).norm() << std::endl;
        std::cout << "Debug end" << std::endl;
#endif

        auto tpl = Get_C_0_Majs(Energys(eval, tr), evec);
        C_0 = std::get<0>(tpl);


        Maj1 = std::get<1>(tpl), Maj2 = std::get<2>(tpl);
    }
    mkl_set_dynamic(0);
    mkl_set_num_threads(MKL_NUM_THREADS);
    omp_set_nested(1);
#ifdef DEBUG_ACTIVE
#pragma omp critical
    {
        std::cout << "Treads = " << OMP_NUM_THREADS << std::endl;
    }
#endif
#pragma omp parallel for shared(Rho_t, C_0, eval, evec, Omegas, std::cout, Maj1, Maj2) default(none) num_threads(OMP_NUM_THREADS)
    for (size_t k = 0; k < Omegas.size(); k++) {
#ifdef DEBUG_ACTIVE
#pragma omp critical
        {
            std::cout << "Tread ID = " << omp_get_thread_num() << std::endl;
        }
#endif
        Vec_cd C_f(MATRIX_SIZE);
        boost::numeric::odeint::bulirsch_stoer<Vec_cd> state;
        boost::numeric::odeint::integrate_const(
                state,
                Schroedinger_of_cs<HAMILTONIAN>(Omegas[k]),
                C_0,
                (double) T_START,

                (double) T_START + T_END,
                double(2.00) * M_PI / T_RES,
                last_observer(C_f));
        double norm = C_f.norm();
        double res = std::pow(std::abs(C_f.dot(Maj1 - Maj2)), 2) / norm / 2;
        double res2 = std::pow(std::abs(C_f.dot(Maj1 + Maj2)), 2) / norm / 2;
        Rho_t.col(k) = Eigen::Vector3d({res, res2, norm});
    }
    return Rho_t;
}
