// Created by Jakob Teuffel on 16.12.19.


#pragma once


#include <complex>
#include <Eigen/Dense>
#include "Hamiltonian_Matrix.hpp"
#include "Typedefs.hpp"

/// This Class generates matrices representing the topological chiral 2D p-Wave Superconducter Hamiltonian
class _2DTopSuperConMatrix : public Hamiltonian_Matrix {
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
    /// Phi
    double phi;
    ///  Helper Funtion to set \f$ A , B |f$ and \f$ \epsilon B^* ,  \f$ simulataneously.
    /// \param k  k - Index
    /// \param j j - Index
    /// \param val New Value of the Submatrix
    inline void set(size_t k, size_t j, cd val) {
        auto MS = 2 * Msize;
        m(k, j) = val;
        m((k + Msize) % MS, (j + Msize) % MS) = -std::conj(val);
    }

    /// Helper function to build the \f$ A \f$ submatrices
    void Build_A() {

<<<<<<< HEAD
        double r = 4*pow(4 / (Gsize_d * 0.05 * phi), 2);
=======
        double r = 4 * pow(4 / (Gsize_d * 0.05), 2);
#pragma unrollandfuse
>>>>>>> 681774e4aadd278002bf4b37513dab9bf4d6039b
        for (size_t k = 0; k < Msize; k++) {
            for (size_t j = 0; j < Msize; j++) {
                size_t x_j = index_x(j);
                size_t y_j = index_y(j);
                size_t x_k = index_x(k);
                size_t y_k = index_y(k);

                if ((x_k == x_j) && (y_k == y_j)) {
                    double Ex = exp(-r * (pow(x_j - Gsize_d / 2 - Vort_x, 2)
                                          + pow(y_j - Gsize_d / 2 - Vort_y, 2)))
                                + exp(-r * (pow(x_j - Gsize_d / 2 + Vort_x, 2)
                                            + pow(y_j - Gsize_d / 2 + Vort_y, 2)));
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
                } else {
                    set(k, j, 0);
                    set(k, j + Msize, 0);
                }
            }
        }
    }

    /// Helper function to build the \f$ B \f$ submatrices
    void Build_B() {
#pragma unrollandfuse
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
                } else {
                    set(k, j + Msize, 0);
                }
            }
        }
    }

    /// Helper fuction to get the index from a x/y coordinate Pair
    /// \param index_x  x - Coordinate
    /// \param index_y  y - Coordinate
    /// \return Index
    inline size_t at(size_t index_x, size_t index_y) {
        return index_x + Gsize * index_y;
    }

    /// Helper Function to get the X - Coordinate from a Index
    /// \param pos - Index
    /// \return X - Coordinate
    inline size_t index_x(size_t pos) {
        return pos % Gsize;
    }

    /// Helper Function to get the Y - Coordinate from a Index
    /// \param pos - Index
    /// \return Y - Coordinate
    inline size_t index_y(size_t pos) {
        return pos / Gsize;
    }

public:
    /// Constructor to create Matrices
    /// \param size_in  Matrix Size
    /// \param t_in tunneling energy
    /// \param phi_in Rotation Angel of Vortice pair
    /// \param mu_in Chemical Potential Amplitude
    /// \param Delta_in Pairing potential
    _2DTopSuperConMatrix(size_t size_in,
                         double t_in,
                         double phi_in,
                         double mu_in,
                         double Delta_in) {
        this->t = t_in;
        this->mu = mu_in;
        this->Delta = static_cast<cd>(Delta_in / 2);
        this->Msize = size_in / 2;
        this->Gsize = sqrt(Msize);
        this->Gsize_d = static_cast<double>(Gsize);
        this->phi =  phi_in;
        {
<<<<<<< HEAD
            //double sp = sin(phi_in), cp = cos(phi_in);
            this->Vort_x = Gsize/4; //sp * (Gsize_d - 1) / 4;
            this->Vort_y = 0.00; //cp * (Gsize_d - 1) / 4;
=======
            double sp = sin(phi_in), cp = cos(phi_in);
            this->Vort_x = sp * (Gsize_d - 1) / 4;
            this->Vort_y = cp * (Gsize_d - 1) / 4;
>>>>>>> 681774e4aadd278002bf4b37513dab9bf4d6039b
        }

#ifdef DEBUG_ACTIVE
#pragma omp critical
        std::cout << "phi= " << phi_in << " Vort_x= " << Vort_x << "  Vort_y = " << Vort_y << std::endl;
#endif
        m = Mat_cd(size_in, size_in);
        Build_A();
        Build_B();
    }

    /// Nessesary override to provide access to the Eigen Storage Matrix
    /// \return Eigen Matrix with coeffitents of Hamiltonian
    Mat_cd get() override {
        return m;
    }

    /// Nessesary override to provide access to the trace
    /// \return trace
    double trace_A() override {
        return 0;
    }
};

