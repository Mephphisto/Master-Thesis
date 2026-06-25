//
// Created by Jakob Teuffel on 08.04.20.
//
#pragma once

#ifdef ITT_ACTIVE

#include <ittnotify.h>

#endif // ITT_ACTIVE

#include <Eigen/Dense>
#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen.hpp>

// Boost.Odeint's vector_space_norm_inf returns complex<T> for complex Eigen
// vectors, but the stepper needs a real scalar. Specialize to fix.
namespace boost { namespace numeric { namespace odeint {
template<typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct vector_space_norm_inf<Eigen::Matrix<std::complex<Scalar>, Rows, Cols, Options, MaxRows, MaxCols>> {
    typedef Scalar result_type;
    result_type operator()(const Eigen::Matrix<std::complex<Scalar>, Rows, Cols, Options, MaxRows, MaxCols>& x) const {
        return x.cwiseAbs().maxCoeff();
    }
};
}}}

#include "Hamiltonian_Matrix.hpp"
#include "Parameters.hpp"
#include "Majorana_Transform.hpp"
#include "Typedefs.hpp"
#include <cmath>
#include <omp.h>
#include <string>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <string>

// #include "_2DTopSuperConMatrix_DEPRECATED.hpp"

typedef double Time;

#if defined(USE_MAGMA) && defined(USE_GPU)
#include "Sim_Magma_Solver.hpp"
typedef magma::SelfAdjointEigenSolver<Mat_cd> MSolver;
#else
typedef Eigen::SelfAdjointEigenSolver<Mat_cd> MSolver;
#endif

inline std::tuple<Vec_cd, Vec_cd, Vec_cd, Vec_cd>
Get_C_0_Majs(const Vec_cd &eval, const Mat_cd &evec, double tr)
{
    // Summ up over Fermmie see and Majorana states
    HAMILTONIAN M(MATRIX_SIZE, T_COUPLE, T_START, MU, DELTA);
    size_t majoranas[2] = {0, 0};
    // Find Majoranas
#pragma unroll
    for (size_t t = 0; t < eval.size(); t++)
    {
        if (abs(eval[t]) < abs(eval[majoranas[0]]))
        {
            majoranas[1] = majoranas[0];
            majoranas[0] = t;
        }
        else if (abs(eval[t]) < abs(eval[majoranas[1]]))
        {
            majoranas[1] = t;
        }
    }
    auto tpl = mjize::FermiiMajize(evec.col(majoranas[0]), evec.col(majoranas[1]));
    return std::make_tuple(std::get<0>(tpl), std::get<1>(tpl), std::get<2>(tpl), std::get<3>(tpl));
}

inline Vec Energys(const Vec_cd &eval, const double &tr)
{
    Vec res = eval.col(0).real() - Vec::Constant(MATRIX_SIZE, 0.5 * (tr - eval.col(0).real().sum()));
    return res;
}

struct Schroedinger_of_cs
{
private:
    double w;
    HAMILTONIAN *H;

public:
    explicit Schroedinger_of_cs(const double &omega, HAMILTONIAN *H_in) : w(omega)
    {
        H = H_in;
#ifdef DEBUG_ACTIVE2
        _2DTopSuperConMatrix_DEPR Mcheck(MATRIX_SIZE, T_COUPLE, T_START, MU, DELTA);
#pragma omp critical
        std::cout << "Diff old new Mat" << (H->get() - Mcheck.get()).norm() << std::endl;
#endif
    }

    /// This is the ODE to be solved
    inline void operator()(const Vec_cd &c, Vec_cd &dcdt, const Time theta)
    {
        H->Update(MATRIX_SIZE, T_COUPLE, theta, MU, DELTA);
        /*Mat_cd M = _2DTopSuperConMatrix_DEPR(MATRIX_SIZE, T_COUPLE, theta, MU, DELTA).get();
        if (std::cout << "Diff Vec " << ((H->get()*c).eval() - (M * c).eval()).norm() << std::numeric_limits<double>::epsilon()){
            std::cout << std::endl << "tid " << omp_get_thread_num << "\t theta " << theta << " \t |c| " << c.norm() << std::endl;
        }
        dcdt = H->getH() * c;*/
        dcdt = H->get() * c; //+ E_offset;//+ H->trace_A() * c;
        // dcdt = M * c;
        dcdt *= cd(0, -1 / w);
    }
};

struct last_observer
{
    Vec_cd &m_state;
    Vec &Progress;
    Eigen::VectorXi &Threads;
    int Threads_old, tid;
    double offset;

    explicit last_observer(Vec_cd &state, Vec &Prog, Eigen::VectorXi &Thread, double offset_in) : m_state(state), Threads(Thread), Progress(Prog), offset(offset_in)
    {
        tid = omp_get_thread_num();
        int Tr_aux;
#pragma omp atomic read
        Tr_aux = Threads[tid];
        Threads_old = Tr_aux;
#ifndef __APPLE__
        mkl_set_num_threads_local(Tr_aux);
#endif
    }

    void operator()(const Vec_cd &x, const double &t)
    {
        m_state = x;

#pragma omp single nowait
        {
            std::cout << " Running " << Progress.size() << " Progress " << Progress.transpose() << std::endl;
        }
        {
#pragma omp atomic write
            Progress[tid] = t / T_END + offset;
            int Tr_aux;
#pragma omp atomic read
            Tr_aux = Threads[tid];
            if (Tr_aux != Threads_old)
            {
#ifndef __APPLE__
                mkl_set_num_threads_local(Tr_aux);
#endif
                Threads_old = Tr_aux;
            }
        }
    }
};

inline void Write_Output_atomic(Mat_cd &Rho_t, size_t k, const Eigen::VectorXcd &res)
{
    for (size_t j = 0; j < res.size(); j++)
    {
#pragma omp atomic write
        reinterpret_cast<double (&)[2]>(Rho_t(j, k))[0] = res[j].real();
#pragma omp atomic write
        reinterpret_cast<double (&)[2]>(Rho_t(j, k))[1] = res[j].imag();
    }
}

inline Eigen::VectorXcd
Prepare_Output(const double Omega, const Vec_cd &Ferm1, const Vec_cd &Ferm2, const Vec_cd &Maj1, const Vec_cd &Maj2,
               double norm, size_t k, const Vec_cd &C_f, double Eergy)
{
    Eigen::VectorXcd res(6);
    res << Omega,
        Ferm1.dot(C_f) / norm,
        Ferm2.dot(C_f) / norm,
        cd(C_f.norm(), Eergy),
        Maj1.dot(C_f) / norm,
        Maj2.dot(C_f) / norm;
    return res;
}

inline void Integrate_Schroedinger(const Vec &Omegas, Mat_cd &Rho_t, const Vec_cd &Ferm1, const Vec_cd &Ferm2,
                                   const Vec_cd &Maj1, const Vec_cd &Maj2, double E_offset, double norm, int tid,
                                   size_t k,
                                   _2DTopSuperConMatrix &myM, Vec_cd &myC_0,
                                   Eigen::VectorXi &Threads,
                                   Vec &Progress, Vec_cd &C_f)
{
    boost::numeric::odeint::bulirsch_stoer<Vec_cd> state;
    // boost::numeric::odeint::runge_kutta_fehlberg78<Vec_cd> state;
    boost::numeric::odeint::integrate_const(
        // boost::numeric::odeint::make_controlled(1e-3,1e-3, state),
        state,
        Schroedinger_of_cs(Omegas[k], &myM),
        myC_0,
        static_cast<double>(T_START),
        static_cast<double>(T_END),
        (T_END - T_START) / T_RES, //(* std::min( Omegas[k]/ E_offset, 1.0),
        last_observer(C_f, Progress, Threads,
                      static_cast<double>(
                          (k < Omegas.size() / 2 ? 2 * (k - tid) : 2 * (Omegas.size() - 1 - k - tid) + 1) *
                          OMP_NUM_THREADS / Omegas.size())));
    {
        myM.Update(MATRIX_SIZE, T_COUPLE, T_END, MU, DELTA);
        double Energy = std::abs(C_f.dot(myM.get() * C_f));
        Write_Output_atomic(Rho_t, k, Prepare_Output(Omegas[k], Ferm1, Ferm2, Maj1, Maj2, norm, k, C_f, Energy));
    }
}

Mat_cd Do_TE(Vec const &Omegas)
{
#ifdef ITT
    __itt_pause();
#endif // ITT
    static_assert(std::is_base_of<Hamiltonian_Matrix, HAMILTONIAN>::value,
                  "Given Class is not derived from Hamiltonian_Matrix");
#ifdef DEBUG_ACTIVE
    std::cout << "MATRIX_SIZE= " << MATRIX_SIZE << std::endl;
#endif
    Mat_cd Rho_t(6, Omegas.size());
    Vec_cd C_0, eval, Ferm1, Ferm2, Maj1, Maj2;
    Mat_cd evec;
    double E_offset, norm;

    HAMILTONIAN M(MATRIX_SIZE, T_COUPLE, T_START, MU, DELTA);

    {
#ifdef DEBUG_ACTIVE
        std::cout << "Hermiticity Check" << std::endl;
        M.verify_hermitiity();
        std::cout << "Matrix Det " << M.get().determinant() << std::endl;
        std::cout << "M*1= " << (M.get().selfadjointView<Eigen::Lower>() * Vec_cd::Ones(MATRIX_SIZE)).norm()
                  << std::endl;
#endif
        double tr = M.trace_A();
        MSolver Solver(MATRIX_SIZE);
        Solver.compute(M.get());
        assert(("Diagonalilisation :: No Success ", Solver.info() == Eigen::Success));
        eval = Solver.eigenvalues();
        evec = Solver.eigenvectors();
        E_offset = eval.cwiseAbs().maxCoeff() / 4;
#ifdef DEBUG_ACTIVE
        std::cout << " Debug Test Diagonalizeing Well?" << std::endl;
        Mat_cd D = (evec.adjoint() * M.get().selfadjointView<Eigen::Lower>() * evec);
        D -= eval.asDiagonal();
        std::cout << "|THT^ -D| = " << D.norm() << std::endl;
        std::cout << "|(THT^ -D)*1| = " << (D.selfadjointView<Eigen::Lower>() * Vec_cd::Ones(MATRIX_SIZE)).norm()
                  << std::endl;
        std::cout << "Debug end" << std::endl;

#endif

        auto tpl = Get_C_0_Majs(eval, evec, tr);

        Maj1 = std::get<0>(tpl).normalized();
        Maj2 = std::get<1>(tpl).normalized();
        Ferm1 = std::get<2>(tpl).normalized();
        Ferm2 = std::get<3>(tpl).normalized();
        C_0 = Ferm1;
        norm = C_0.norm();
#ifdef DEBUG_ACTIVE
        std::cout << "< Fermi1|C_0>" << std::norm(Ferm1.dot(C_0)) / norm << std::endl
                  << "< Fermi2|C_0>" << std::norm(Ferm2.dot(C_0)) / norm << std::endl
                  << "< C_0|C_0>" << C_0.norm() << std::endl
                  << "< Maj1|C_0>" << std::norm(Maj1.dot(C_0)) / norm << std::endl
                  << "< Maj2|C_0>" << std::norm(Maj2.dot(C_0)) / norm << std::endl
                  << "< Fermi1|Fermi2>" << std::norm(Ferm1.dot(Ferm2)) / norm << std::endl
                  << "< Maj1|Maj2>" << std::norm(Maj1.dot(Maj2)) / norm << std::endl;

#endif
#ifdef DEBUG_ACTIVE
        std::cout << "<C_0|M|C_0> = " << C_0.dot(M.get() * C_0) << std::endl;
#endif
    }
#ifndef __APPLE__
    mkl_set_dynamic(0);
#endif
    omp_set_max_active_levels(4);
#ifdef DEBUG_ACTIVE
#pragma omp critical
    {
        std::cout << "Treads = " << OMP_NUM_THREADS << std::endl;
    }
#endif
#ifdef ITT_ACTIVE
    __itt_resume();
#endif // ITT_ACTIVE
    Eigen::VectorXi Threads = (MKL_NUM_THREADS / OMP_NUM_THREADS) * Eigen::VectorXi::Ones(OMP_NUM_THREADS);
    Vec Progress = Vec::Zero(OMP_NUM_THREADS);
#ifndef __APPLE__
    mkl_set_dynamic(0);
    mkl_set_dynamic(0);
#endif
    omp_set_max_active_levels(2);
#pragma omp parallel shared(Rho_t, Threads, Progress, C_0, M, std::cout) firstprivate(Omegas, eval, evec, Ferm1, Ferm2, Maj1, Maj2, norm, E_offset) default(none) num_threads(OMP_NUM_THREADS)
    {
        int tid = omp_get_thread_num();
        for (size_t k = tid; k < Omegas.size() / 2; k += OMP_NUM_THREADS)
        {

#ifdef DEBUG_ACTIVE
#pragma omp critical
            std::cout << "HL_Thread " << omp_get_thread_num() << " of " << omp_get_num_threads() << std::endl;
#endif
            HAMILTONIAN myM(M);
            Vec_cd myC_0 = C_0;
            Vec_cd C_f(MATRIX_SIZE);
            Integrate_Schroedinger(Omegas, Rho_t, Ferm1, Ferm2, Maj1, Maj2, E_offset, norm, tid, k, myM, myC_0, Threads,
                                   Progress, C_f);
            myM = M;
            myC_0 = C_0;
            // C_f(MATRIX_SIZE);
            Integrate_Schroedinger(Omegas, Rho_t, Ferm1, Ferm2, Maj1, Maj2, E_offset, norm, tid, Omegas.size() - k - 1,
                                   myM, myC_0, Threads, Progress, C_f);
        }
#ifndef __APPLE__
        mkl_set_num_threads_local(0);
#endif
        {
            int tr;
#pragma atomic read
            tr = Threads[tid];
#pragma atomic write
            Threads[tid] = 0;
            {
                size_t j;
#pragma atomic
                Threads.minCoeff(&j);
#pragma atomic write
                Threads[j] += tr;
#pragma atomic write
                Threads[tid] = std::numeric_limits<int>::infinity();
            }
        }
    }

#ifdef ITT_ACTIVE
    __itt_detach();
#endif // ITT_ACTIVE
    return Rho_t;
}
