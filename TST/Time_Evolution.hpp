//
// Created by Jakob Teuffel on 08.04.20.
//
#pragma once

#ifdef ITT_ACTIVE

#include <ittnotify.h>

#endif //ITT_ACTIVE

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
#include <string>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>
#include <string>

#include "_2DTopSuperConMatrix_DEPRECATED.hpp"

typedef double Time;


#if  defined(USE_MAGMA) && defined(USE_GPU)
#include "TST_Magma_Solver.hpp"
typedef magma::SelfAdjointEigenSolver<Mat_cd> MSolver;
#else
typedef Eigen::SelfAdjointEigenSolver<Mat_cd> MSolver;
#endif

inline std::tuple<Vec_cd, Vec_cd, Vec_cd> Get_C_0_Majs(const Vec_cd &eval, const Mat_cd &evec, double tr) {
    Vec_cd C_0 = Vec_cd(MATRIX_SIZE);
    // Summ up over Fermmie see and Majorana states
    HAMILTONIAN M(MATRIX_SIZE, T_COUPLE, T_START, MU, DELTA);
    size_t majoranas[2] = {0, 0};
    //Find Majoranas
#pragma unroll
    for (size_t t = 0; t < eval.size(); t++) {
        if (abs(eval[t]) < abs(eval[majoranas[0]])) {
            majoranas[1] = majoranas[0];
            majoranas[0] = t;
        } else if (abs(eval[t]) < abs(eval[majoranas[1]])) {
            majoranas[1] = t;
        }
    }
    /*
#ifdef DEBUG_ACTIVE
    double E = 0;
#endif
#pragma unroll
    for (size_t k = 0; k < eval.size(); k++) {
        if ((k != majoranas[0]) || k > (k != majoranas[1])) {
            const double e_k = eval[k].real();
            if (std::signbit(e_k)) {
                C_0 += evec.col(k);
#ifdef DEBUG_ACTIVE
                E += e_k;
#endif
            }
+
    }
#ifdef DEBUG_ACTIVE
    std::cout << "E_maj1 = " << eval[majoranas[0]] << " E_maj2 = " << eval[majoranas[1]] << std::endl;
#endif
    auto tpl = Majoranaize(evec.col(majoranas[0]).normalized(), evec.col(majoranas[1]).normalized());
#ifdef DEBUG_ACTIVE
    std::cout << "|C_0| = " << C_0.norm() << std::endl;
#endif
    C_0 = (C_0 + ((0 < eval[majoranas[0]].real()) ? evec.col(majoranas[1]) : evec.col(majoranas[0]))).eval();
    //C_0 += evec.col(majoranas[0]).normalized() + evec.col(majoranas[1]).normalized();
    //C_0 += std::get<0>(tpl) ; //+ std::get<1>(tpl);
#ifdef DEBUG_ACTIVE
    E += std::max(eval[majoranas[0]].real(), eval[majoranas[1]].real());
    std::cout << "E_0 = " << E << std::endl;
    double a = std::abs(E - (C_0.dot(M.get() * C_0) / C_0.squaredNorm()));
    if (a > 1e-12) std::cout << " E-Hv/vv = " << a << " E= "<<E << " Hv/vv= " << (C_0.dot(M.get() * C_0)) <<  std::endl;
#endif
    return std::make_tuple(C_0, std::get<0>(tpl), std::get<1>(tpl));*/
    C_0 = evec.col(majoranas[0]);
    return std::make_tuple(C_0, evec.col(majoranas[0]), evec.col(majoranas[1]));
}

inline Vec Energys(const Vec_cd &eval, const double &tr) {
    Vec res = eval.col(0).real() - Vec::Constant(MATRIX_SIZE, 0.5 * (tr - eval.col(0).real().sum()));
    return res;
}


struct Schroedinger_of_cs {
private:
    double w, E_offset;
    HAMILTONIAN *H;
public:
    explicit Schroedinger_of_cs(const double &omega, const double &E_ofs, HAMILTONIAN *H_in) : w(omega),
                                                                                               E_offset(E_ofs) {
        H = H_in;
#ifdef DEBUG_ACTIVE2
        _2DTopSuperConMatrix_DEPR Mcheck(MATRIX_SIZE, T_COUPLE, T_START, MU, DELTA);
#pragma omp critical
        std::cout << "Diff old new Mat" << (H->get() - Mcheck.get()).norm() << std::endl;
#endif
    }

    /// This is the ODE to be solved
    inline void operator()(const Vec_cd &c, Vec_cd &dcdt, const Time theta) {
        H->Update(MATRIX_SIZE, T_COUPLE, theta, MU, DELTA);
        /*Mat_cd M = _2DTopSuperConMatrix_DEPR(MATRIX_SIZE, T_COUPLE, theta, MU, DELTA).get();
        if (std::cout << "Diff Vec " << ((H->get()*c).eval() - (M * c).eval()).norm() << std::numeric_limits<double>::epsilon()){
            std::cout << std::endl << "tid " << omp_get_thread_num << "\t theta " << theta << " \t |c| " << c.norm() << std::endl;
        }
        dcdt = H->getH() * c;*/
        dcdt = H->get() * c + E_offset;//+ H->trace_A() * c;
        //dcdt = M * c;
        dcdt *= cd(0, -1 / w);
    }
};

struct last_observer {
    Vec_cd &m_state;
    Vec &Progress;
    Eigen::VectorXi &Threads;
    int Threads_old, tid;
    double offset;

    explicit last_observer(Vec_cd &state, Vec &Prog, Eigen::VectorXi &Thread, double offset_in) :
            m_state(state), Threads(Thread), Progress(Prog), offset(offset_in) {
        tid = omp_get_thread_num();
        int Tr_aux;
#pragma omp atomic read
        Tr_aux = Threads[tid];
        Threads_old = Tr_aux;
        mkl_set_num_threads_local(Tr_aux);
    }

    void operator()(const Vec_cd &x, const double &t) {
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
            if (Tr_aux != Threads_old) {
                mkl_set_num_threads_local(Tr_aux);
                Threads_old = Tr_aux;
            }
        }
    }
};

Mat Do_TE(Vec const &Omegas) {
#ifdef ITT
    __itt_pause();
#endif //ITT
    static_assert(std::is_base_of<Hamiltonian_Matrix, HAMILTONIAN>::value,
                  "Given Class is not derived from Hamiltonian_Matrix");
#ifdef DEBUG_ACTIVE
    std::cout << "MATRIX_SIZE= " << MATRIX_SIZE << std::endl;
#endif
    Mat Rho_t(6, Omegas.size());
    Vec_cd C_0, eval, Maj1, Maj2;
    Mat_cd evec;
    double E_offset;

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
        E_offset = eval.cwiseAbs().maxCoeff() * (R_STIFF + 1.0) / (R_STIFF - 1.0);
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
        C_0 = std::get<0>(tpl);
#ifdef DEBUG_ACTIVE
        std::cout << "<c_0|H|c_0> = " << C_0.dot(M.get() * C_0) + tr << std::endl;
#endif

        Maj1 = std::get<1>(tpl), Maj2 = std::get<2>(tpl);
    }
    double norm = std::pow(std::abs(C_0.dot(Maj1 + Maj2)), 2);
    mkl_set_dynamic(0);
    omp_set_max_active_levels(4);
#ifdef DEBUG_ACTIVE
#pragma omp critical
    {
        std::cout << "Treads = " << OMP_NUM_THREADS << std::endl;
    }
#endif
#ifdef ITT_ACTIVE
    __itt_resume();
#endif //ITT_ACTIVE
    Eigen::VectorXi Threads = (MKL_NUM_THREADS / OMP_NUM_THREADS) * Eigen::VectorXi::Ones(OMP_NUM_THREADS);
    Vec Progress = Vec::Zero(OMP_NUM_THREADS);
    mkl_set_dynamic(0);
    omp_set_max_active_levels(2);
#pragma omp parallel shared(Rho_t, Threads, Progress, C_0, M, std::cout) firstprivate(Omegas, eval, evec, Maj1, Maj2, norm, E_offset) default(none) num_threads(OMP_NUM_THREADS)
    {
        int tid = omp_get_thread_num();
        for (size_t k = tid; k < Omegas.size() / 2; k += OMP_NUM_THREADS) {

#ifdef DEBUG_ACTIVE
#pragma omp critical
            std::cout << "HL_Thread " << omp_get_thread_num() << " of " << omp_get_num_threads() << std::endl;
#endif
            Vec_cd myC_0 = C_0;
            HAMILTONIAN myM(M);
            Vec_cd C_f(MATRIX_SIZE);
            boost::numeric::odeint::bulirsch_stoer<Vec_cd> state;
            boost::numeric::odeint::integrate_const(
                    state,
                    Schroedinger_of_cs(Omegas[k], E_offset, &myM),
                    myC_0,
                    static_cast<double>(T_START),
                    static_cast<double>(T_END),
                    static_cast<double>(2.00) * M_PI / T_RES / (R_STIFF + 1.0),//* Omegas[k],
                    last_observer(C_f, Progress, Threads, static_cast<double>(k - tid) * 2 / Omegas.size()));
            {
                Eigen::VectorXd res(6);
                res << Omegas[k],
                        std::norm(C_f.dot(Maj1)) / norm,
                        std::norm(C_f.dot(Maj2)) / norm,
                        C_f.norm(),
                        std::arg((C_f.transpose() * Maj1)[0]),
                        std::arg((C_f.transpose() * Maj2)[0]);

                for (size_t j = 0; j < 5; j++) {
#pragma omp atomic write
                    Rho_t(j, k) = res[j];
                }
            }
            myC_0 = C_0;
            myM = M;
            //C_f(MATRIX_SIZE);
            boost::numeric::odeint::integrate_const(
                    state,
                    Schroedinger_of_cs(Omegas[Omegas.size() - k - 1], E_offset, &myM),
                    myC_0,
                    static_cast<double>(T_START),
                    static_cast<double>(T_END),
                    static_cast<double>(2.00) * M_PI / T_RES / (R_STIFF + 1.0),//* Omegas[k],
                    last_observer(C_f, Progress, Threads, (static_cast<double>(k - tid) * 2 + 1) / (Omegas.size())));
            {
                Eigen::VectorXd res(6);
                res << Omegas[Omegas.size() - k - 1],
                        std::pow(std::abs(C_f.dot(Maj1 - Maj2)), 2) / norm,
                        std::pow(std::abs(C_f.dot(Maj1 + Maj2)), 2) / norm,

                        C_f.norm(),
                        std::arg((C_f.transpose() * Maj1)[0]),
                        std::arg((C_f.transpose() * Maj2)[0]);

                for (size_t j = 0; j < 5; j++) {
#pragma omp atomic write
                    Rho_t(j, Omegas.size() - k - 1) = res[j];
                }
            }
        }
        mkl_set_num_threads_local(0);
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
#endif //ITT_ACTIVE
    return
            Rho_t;
}
