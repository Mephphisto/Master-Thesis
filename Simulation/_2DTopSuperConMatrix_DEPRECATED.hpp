// Created by Jakob Teuffel on 16.12.19.


#pragma once


#include <complex>
#include <Eigen/Dense>
#include "Hamiltonian_Matrix.hpp"
#include "Typedefs.hpp"

/// DEPRECATED Version of _2DTopSuperConMatrix cn be used to verify Optimized / obfuctaded new code


class _2DTopSuperConMatrix_DEPR : public Hamiltonian_Matrix {
private:
    /// Tunneling Energy
    double t;
    /// Chemical Potential Amplitude
    double mu;
    /// Grid Size as double for computations
    double Gsize_d;
    /// Vortex X-Position
    double Vort_x;
    /// Vortex Y-Position
    double Vort_y;
    /// Complex Pairing Potential
    cd Delta;
    /// Grid  size
    size_t Gsize;

    ///  Helper Funtion to set \f$ A , B |f$ and \f$ \epsilon B^* ,  \f$ simulataneously.
    /// \param k  k - Index
    /// \param j j - Index
    /// \param val New Value of the Submatrix
    inline void set(const size_t &k, const size_t &j, const cd &val) {
        size_t MS = 2 * Msize;
        size_t j2 = (j + Msize) % MS;
        size_t k2 = (k + Msize) % MS;
        cd val2 = -std::conj(val);
        m(k, j) = val;
        m(k2, j2) = val2;
    }

    /// Helper function to build the \f$ A \f$ submatrices
    inline void Build_A() {
        double r = 4 * pow(4 / (Gsize_d * 0.05), 2);

#ifdef ICC
#pragma unroll_and_jam
#else
#pragma unrollandfuse
#endif
        for (size_t k = 0; k < Msize; k++) {
            for (size_t j = 0; j < Msize; j++) {
                size_t x_j = index_x(j);
                size_t y_j = index_y(j);
                size_t x_k = index_x(k);
                size_t y_k = index_y(k);

                if (x_k == x_j) {
                    if (y_k == y_j) {
                        double Ex = exp(-r * (pow(x_j - Gsize_d / 2 - Vort_x, 2)
                                              + pow(y_j - Gsize_d / 2 - Vort_y, 2)))
                                    + exp(-r * (pow(x_j - Gsize_d / 2 + Vort_x, 2)
                                                + pow(y_j - Gsize_d / 2 + Vort_y, 2)));
                        set(k, j, mu - 2 * mu * Ex);

                    } else if (y_k == y_j + 1) {
                        set(k, j, t);
                    } else if (y_k == y_j - 1) {
                        set(k, j, t);
                    }
#ifdef PERIODIC_BOUNDRY
                    else if ((y_k == Gsize - 1) && (y_j == 0)) {
                        set(k, j, t);
                    } else if ((y_k == 0) && (y_j == Gsize - 1)) {
                        set(k, j, t);
                    }
#endif
                } else if (y_k == y_j) {
                    if (x_k == x_j + 1) {
                        set(k, j, t);
                    } else if (x_k == x_j - 1) {
                        set(k, j, t);
                    }
                }
#ifdef PERIODIC_BOUNDRY
                if ((x_k == Gsize - 1) && (x_j == 0) && (y_k == y_j)) {
                   set(k, j, t);
               } else if ((x_k == 0) && (x_j == Gsize - 1) && (y_k == y_j)) {
                   set(k, j, t);
               }
#endif
                /*else {
                   set(k, j, 0);
                   set(k, j + Msize, 0);
               }*/
            }
        }
    }


    cd inline getPhase(const size_t &x_j, const size_t &y_j, const size_t &x_k, const size_t &y_k) const {
        cd phase;
#ifdef PHASE

        phase = cd(1, 0);
        phase *= cd(static_cast<double>(x_j + x_k) / 2 - (Gsize_d - 1) / 2 - Vort_x,
                    static_cast<double>(y_j + y_k) / 2 - (Gsize_d - 1) / 2 - Vort_y);
        phase *= cd(static_cast<double>(x_j + x_k) / 2 - (Gsize_d - 1) / 2 + Vort_x,
                    static_cast<double>(y_j + y_k) / 2 - (Gsize_d - 1) / 2 + Vort_y);
        {
            double Abs = abs(phase);
            if (Abs != 0.0) {
                phase = phase / Abs;
            } else {
                phase = 0.0;
            }
        }
#else
        phase = cd(1, 0);
#endif
        return phase;
    }

/// Helper function to build the \f$ B \f$ submatrices
    inline void Build_B() {

#ifdef ICC
#pragma unroll_and_jam
#else
#pragma unrollandfuse
#endif
        for (size_t k = 0; k < Msize; k++) {
            for (size_t j = 0; j < Msize; j++) {
                size_t x_j = index_x(j);
                size_t y_j = index_y(j);
                size_t x_k = index_x(k);
                size_t y_k = index_y(k);

                if (x_k == x_j) {
                    if (y_k == y_j + 1) {
                        set(k, j + Msize, cd(0, 1) * Delta * getPhase(x_j, y_j, x_k, y_k));
                    } else if (1 + y_k == y_j) {
                        set(k, j + Msize, cd(0, -1) * Delta * getPhase(x_j, y_j, x_k, y_k));
                    }
#ifdef PERIODIC_BOUNDRY
                    else if ((y_k == Gsize - 1) && (y_j == 0)) {
                        set(k, j + Msize, cd(0, -1) * Delta );
                    } else if ((y_k == 0) && (y_j == Gsize - 1)) {
                        set(k, j + Msize, cd(0, 1) * Delta);
                    }
#endif
                } else if (y_k == y_j) {
                    if (x_k == x_j + 1) {
                        set(k, j + Msize, cd(1, 0) * Delta * getPhase(x_j, y_j, x_k, y_k));
                    } else if (1 + x_k == x_j) {
                        set(k, j + Msize, cd(-1, 0) * Delta * getPhase(x_j, y_j, x_k, y_k));
                    }
                }
#ifdef PERIODIC_BOUNDRY
                if ((x_k == Gsize - 1) && (x_j == 0) && (y_k == y_j)) {
                   set(k, j + Msize, cd(-1, 0) * Delta );
               } else if ((x_k == 0) && (x_j == Gsize - 1) && (y_k == y_j)) {
                   set(k, j + Msize, cd(1, 0) * Delta );
               }
#endif
                /* else {
                set(k, j + Msize, 0);
            }*/
            }
        }
    }


public:
    /// Helper fuction to get the index from a x/y coordinate Pair
    /// \param index_x  x - Coordinate
    /// \param index_y  y - Coordinate
    /// \return Index
    inline size_t at(const size_t &index_x, const size_t &index_y) {
        return index_x + Gsize * index_y;
    }

    /// Helper Function to get the X - Coordinate from a Index
    /// \param pos - Index
    /// \return X - Coordinate
    inline size_t index_x(const size_t &pos) {
        return pos % Gsize;
    }

    /// Helper Function to get the Y - Coordinate from a Index
    /// \param pos - Index
    /// \return Y - Coordinate
    inline size_t index_y(size_t pos) {
        return pos / Gsize;
    }


    /// Constructor to create Matrices
    /// \param size_in  Matrix Size
    /// \param t_in tunneling energy
    /// \param phi_in Rotation Angel of Vortice pair
    /// \param mu_in Chemical Potential Amplitude
    /// \param Delta_in Pairing potential
    _2DTopSuperConMatrix_DEPR(
            const size_t &size_in,
            const double &t_in,
            const double &phi_in,
            const double &mu_in,
            const double &Delta_in) :
            t(t_in),
            mu(mu_in),
            Hamiltonian_Matrix(size_in),
            Delta(static_cast<cd>(Delta_in / 2)) {
        this->Msize = size_in / 2;
        this->Gsize = sqrt(Msize);
        this->Gsize_d = static_cast<double>(Gsize);
        {
            double sp = sin(phi_in), cp = cos(phi_in);
            this->Vort_x = sp * (Gsize_d - 1) / 4;
            this->Vort_y = cp * (Gsize_d - 1) / 4;
        }

#ifdef DEBUG_ACTIVE2
#pragma omp critical
        std::cout << "phi= " << phi_in << " Vort_x= " << Vort_x << "  Vort_y = " << Vort_y << std::endl;
#endif
        m = Mat_cd::Zero(size_in, size_in);
        Build_A();
        Build_B();
    }

    /// Update Matrix
    /// \param size_in  Matrix Size
    /// \param t_in tunneling energy
    /// \param phi_in Rotation Angel of Vortice pair
    /// \param mu_in Chemical Potential Amplitude
    /// \param Delta_in Pairing potential
    void inline Update(
            const size_t &size_in,
            const double &t_in,
            const double &phi_in,
            const double &mu_in,
            const double &Delta_in) {
        this->t = t_in;
        this->mu = mu_in;
        this->Delta = static_cast<cd>(Delta_in / 2);
        this->Msize = size_in / 2;
        this->Gsize = sqrt(Msize);
        this->Gsize_d = static_cast<double>(Gsize);
        {
            double sp = sin(phi_in), cp = cos(phi_in);
            this->Vort_x = sp * (Gsize_d - 1) / 4;
            this->Vort_y = cp * (Gsize_d - 1) / 4;
        }

#ifdef DEBUG_ACTIVE2
#pragma omp critical
        std::cout << "phi= " << phi_in << " Vort_x= " << Vort_x << "  Vort_y = " << Vort_y << std::endl;
#endif
        Build_A();
        Build_B();
    }


    /// Nessesary override to provide access to the trace
    /// \return trace
    double trace_A()
    override {
        return 0;
    }
};

