// Created by Jakob Teuffel on 16.12.19.


#pragma once

#include "Hamiltonian_MatrixSparse.hpp"
#include "Typedefs.hpp"
#include <complex>
#include <omp.h>
#include <Eigen/Sparse>
/// This Class is a deprecated attempt to convert the Prohject to sparse matrices. This Class is however functional
#pragma deprecated
class _2DTopSuperConMatrixSparse : public Hamiltonian_MatrixSparse {
private:
    double t, mu, Gsize_d, Vort_x , Vort_y  ;
    cd Delta;
    size_t Gsize;

    inline void set(size_t k, size_t j, cd val) {
        auto MS = 2 * Msize;
        m.insert(k, j) = val;
        m.insert((k + Msize) % MS, (j + Msize) % MS) = -std::conj(val);
    }

    inline void Build_A() {

        double r = pow(4/(Gsize_d*0.05),2);
#pragma unrollandfuse
        for (size_t k = 0; k < Msize; k++) {
            for (size_t j = 0; j < Msize; j++) {
                size_t x_j = index_x(j);
                size_t y_j = index_y(j);
                size_t x_k = index_x(k);
                size_t y_k = index_y(k);

                if ((x_k == x_j) && (y_k == y_j)) {
                    double Ex = exp(-r * (pow(x_j - Gsize_d / 2 - Vort_x, 2)
                                          +  pow(y_j - Gsize_d / 2 - Vort_y, 2)))
                                + exp(-r * (pow(x_j - Gsize_d / 2 + Vort_x, 2)
                                            +  pow(y_j - Gsize_d / 2 + Vort_y, 2)));
                    set(k, j, mu - 2 * mu * Ex);
                    // Alternative Code with Hard Edge
                    /*if
                    (sqrt(pow(x_j - Gsize_d / 2 - Vort_x, 2) + pow(y_j - Gsize_d / 2 - Vort_y, 2)) < .3 * Gsize_d / 4) {
                        set(k, j, -mu);
                    } else if
                    (sqrt(pow(x_j - Gsize_d / 2 + Vort_x, 2) + pow(y_j - Gsize_d / 2 + Vort_y, 2)) < .3 * Gsize_d / 4) {
                        set(k, j, - mu);
                    } else {
                        set(k, j,mu);
                }*/

                } else if ((x_k == x_j) && (y_k == y_j + 1)) {
                    set(k, j, t);
                } else if ((x_k == x_j + 1) && (y_k == y_j)) {
                    set(k, j, t);
                } else if ((x_k == x_j) && (y_k == y_j - 1)) {
                    set(k, j, t);
                } else if ((x_k == x_j - 1) && (y_k == y_j)) {
                    set(k, j, t);
#ifdef PERIODIC_BOUNDRY
                    } else if ((x_k == x_j) && (y_k == Gsize - 1) && (y_j == 0)) {
                        set(k, j, t);
                    } else if ((x_k == x_j) && (y_k == 0) && (y_j == Gsize - 1)) {
                        set(k, j, t);
                    } else if ((x_k == Gsize - 1) && (x_j == 0) && (y_k == y_j)) {
                        set(k, j, t);
                    } else if ((x_k == 0) && (x_j == Gsize - 1) && (y_k == y_j)) {
                        set(k, j, t);
#endif
                }
            }
        }
    }

    inline void Build_B() {
#pragma unroll
        for (size_t k = 0; k < Msize; k++) {
            for (size_t j = 0; j < Msize; j++) {
                size_t x_j = index_x(j);
                size_t y_j = index_y(j);
                size_t x_k = index_x(k);
                size_t y_k = index_y(k);

                cd phase;
#ifdef PHASE

                phase = cd(1, 0);
                phase *= cd(static_cast<double>(x_j + x_k) / 2 - (Gsize_d - 1) / 2 - Vort_x,
                                              static_cast<double>(y_j + y_k) / 2 - (Gsize_d - 1) / 2 - Vort_y);
                phase *= cd(static_cast<double>(x_j + x_k) / 2 - (Gsize_d - 1) / 2 + Vort_x,
                                              static_cast<double>(y_j + y_k) / 2 - (Gsize_d - 1) / 2 + Vort_y);
                phase = phase / abs(phase);
#else
                phase = cd(1,0);
#endif
                if ((x_k == x_j) && (y_k == y_j + 1)) {
                    set(k, j + Msize, cd(0, 1) * Delta * phase);
                } else if ((x_k == x_j) && (y_k == y_j - 1)) {
                    set(k, j + Msize, cd(0, -1) * Delta * phase);
                } else if ((x_k == x_j + 1) && (y_k == y_j)) {
                    set(k, j + Msize, cd(1, 0) * Delta * phase);
                } else if ((x_k == x_j - 1) && (y_k == y_j)) {
                    set(k, j + Msize, cd(-1, 0) * Delta * phase);
#ifdef PERIODIC_BOUNDRY
                    } else if ((x_k == x_j) && (y_k == Gsize - 1) && (y_j == 0)) {
                        set(k, j + Msize, cd(0, -1) * Delta * phase);
                    } else if ((x_k == x_j) && (y_k == 0) && (y_j == Gsize - 1)) {
                        set(k, j + Msize, cd(0, 1) * Delta * phase);
                    } else if ((x_k == Gsize - 1) && (x_j == 0) && (y_k == y_j)) {
                        set(k, j + Msize, cd(-1, 0) * Delta * phase);
                    } else if ((x_k == 0) && (x_j == Gsize - 1) && (y_k == y_j)) {
                        set(k, j + Msize, cd(1, 0) * Delta * phase);
#endif
                }
            }
        }
    }

    inline size_t at(const size_t &index_x, const size_t &index_y) {
        return index_x + Gsize * index_y;
    }

    inline size_t index_x(const size_t &pos) {
        return pos % Gsize;
    }

    inline size_t index_y(const size_t &pos) {
        return pos / Gsize;
    }

public:
    _2DTopSuperConMatrixSparse(const size_t &size_in,
                               const double &t_in,
                               const double &phi_in,
                               const double &mu_in,
                               const double &Delta_in) :
            t(t_in),
            mu(mu_in),
            Delta(static_cast<cd>(Delta_in / 2)),
            Gsize(sqrt(Msize)),
            Gsize_d(static_cast<double>(Gsize)),
            Vort_x(sin(phi_in) * (Gsize_d - 1) / 4),
            Vort_y(cos(phi_in) * (Gsize_d - 1) / 4) {
        this->Msize = size_in / 2;
        {
            // send message in a unified Racecondition safe way
            std::string thread_msg = (omp_get_num_threads() > 1) ? "of thread " + (std::to_string(omp_get_thread_num()))
                                                                 : "";
            std::string msg = "Vortex" + thread_msg +
                              " at (" + std::to_string(Vort_x) + " , " + std::to_string(Vort_y) + " ) \n";
            std::cout << msg;
        }
        m.resize(size_in, size_in);
        m.reserve(size_in * 9);
        Build_A();
        Build_B();
    }

    Eigen::SparseMatrix<cd> &get() override {
        return m;
    }

    double trace_A() override {
        return 0;
    }
};

