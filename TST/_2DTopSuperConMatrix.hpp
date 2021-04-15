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

    ///  Helper Funtion to set \f$ A , B |f$ and \f$ \epsilon B^* ,  \f$ simulataneously.
    /// \param k  k - Index
    /// \param j j - Index
    /// \param val New Value of the Submatrix
//#pragma declare simd
    inline void set(const size_t &k, const size_t &j, const cd &val) {
        m(k, j) = val;
        const size_t Msize2 = 2 * Msize;
        m((k + Msize) % Msize2, (j + Msize) % Msize2) = -std::conj(val);
    }

    /// Helper function to build the \f$ A \f$ submatrices
    inline void Build_A() {
         double r = -std::log(0.5) / RADIUS / RADIUS;

#ifdef ICC
#pragma unroll_and_jam
#else
#pragma unrollandfuse
#endif
        for (size_t x = 0; x < Gsize; x++) {
            double E1x = exp(-r * (pow(x - Gsize_d / 2 - Vort_x, 2))),
                    E2x = exp(-r * (pow(x - Gsize_d / 2 + Vort_x, 2)));
//#pragma omp simd
            for (size_t y = 1; y < Gsize; y++) {
                //double Ex = E1x * exp(-r * pow(y - Gsize_d / 2 - Vort_y, 2))
                //            + E2x * exp(-r * pow(y - Gsize_d / 2 + Vort_y, 2));
                //set(at(x, y), at(x, y), mu - 2 * mu * Ex);
                set(at(x, y), at(x, y),
                    -4 * t + (r < pow(x - Gsize_d / 2 - Vort_x, 2) + pow(y - Gsize_d / 2 + Vort_y, 2)) ? mu : -mu);

                set(at(x, y - 1), at(x, y), t);
                set(at(x, y), at(x, y - 1), t);
                set(at(y, x), at(y - 1, x), t);
                set(at(y - 1, x), at(y, x), t);
            }
            double Ex = E1x * exp(-r * pow(-Gsize_d / 2 - Vort_y, 2))
                        + E2x * exp(-r * pow(-Gsize_d / 2 + Vort_y, 2));
            //set(at(x, 0), at(x, 0), mu - 2 * mu * Ex);
            set(at(x, 0), at(x, 0),
                -4 * t + (r < pow(x - Gsize_d / 2 - Vort_x, 2) + pow(-Gsize_d / 2 + Vort_y, 2)) ? mu : -mu);
        }

    }

//#pragma omp declare simd
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
#pragma unrollandfuse
#endif
        for (size_t x = 0; x < Gsize; x++) {
//#pragma omp simd
            for (size_t y = 0; y < Gsize - 1; y++) {

                set(at(x, y + 1), at(x, y) + Msize, cd(0, 1) * Delta * getPhase(x, y, x, y + 1));
                set(at(x, y), at(x, y + 1) + Msize, cd(0, -1) * Delta * getPhase(x, y + 1, x, y));
                set(at(y + 1, x), at(y, x) + Msize, cd(1, 0) * Delta * getPhase(y, x, y + 1, x));
                set(at(y, x), at(y + 1, x) + Msize, cd(-1, 0) * Delta * getPhase(y + 1, x, y, x));
            }
        }
    }

public:
    /// Helper fuction to get the index from a x/y coordinate Pair
    /// \param index_x  x - Coordinate
    /// \param index_y  y - Coordinate
    /// \return Index
//
//
//
//
//
//
// #pragma omp declare simd
    inline size_t at(const size_t &index_x, const size_t &index_y) {
        return index_x + Gsize * index_y;
    }

    /// Helper Function to get the X - Coordinate from a Index
    /// \param po_NEWs - Index
    /// \return X - Coordinate
    inline size_t index_x(const size_t &pos) const {
        return pos % Gsize;
    }

    /// Helper Function to get the Y - Coordinate from a Index
    /// \param pos - Index
    /// \return Y - Coordinate
    inline size_t index_y(size_t pos) const {
        return pos / Gsize;
    }


    /// Constructor to create Matrices
    /// \param size_in  Matrix Size
    /// \param t_in tunneling energy
    /// \param phi_in Rotation Angel of Vortice pair
    /// \param mu_in Chemical Potential Amplitude
    /// \param Delta_in Pairing potential
    _2DTopSuperConMatrix(const size_t &size_in, const double &t_in, const double &phi_in, const double &mu_in,
                         const double &Delta_in) :
            Hamiltonian_Matrix(size_in),
            t(t_in),
            mu(mu_in),
            Delta(static_cast<cd>(Delta_in / 2)) {
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

    _2DTopSuperConMatrix(const _2DTopSuperConMatrix &H) : Hamiltonian_Matrix(H),
                                                          t(H.t), mu(H.mu), Delta(H.Delta) {
#ifdef DEBUG_ACTIVE
#pragma omp critical
        std::cout << " Copy _2DTopSuperConMatrix " << std::endl;
#endif
        this->Gsize = H.Gsize;
        this->Gsize_d = H.Gsize_d;
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
        return std::abs(m.topLeftCorner<MATRIX_SIZE / 2, MATRIX_SIZE / 2>().trace());
    }
};