


#pragma once

#include "Hamiltonian_Matrix.hpp"
#include "Typedefs.hpp"
#include <Eigen/Dense>

/// This Class generates matrices representing the topological 1D Kitaev Hamiltonian
class _1DTopSuperConMatrix : Hamiltonian_Matrix {
public:
    /// Constructor for the Hamiltonian
    /// \param size_in Matrix size
    /// \param t_con Tunneling Energy
    /// \param param Parameter to be sweapt
    /// \param mu_in Chemical Potential
    /// \param Delta_in Pairing Potential
    _1DTopSuperConMatrix(size_t size_in,
                         double t_con,
                         double param,
                         double mu_in,
                         double Delta_in) :
            Hamiltonian_Matrix(size_in),
            t(-t_con),
            mu(-mu_in),
            Delta(Delta_in / 2),
            param(param),
            Gsize(size_in / 2) {
        Build_A();
    }

    /// Nessesary override to provide access to the trace
    /// \return trace
    inline double trace_A() override {
        return 0;
    }

private:
    /// test a , b , c
    double t, mu, Delta, param;
    /// Gsize is half the Matrix size and thus the size of the \f$ A, B\f$ submatrices
    size_t Gsize;

    /// Sigmoid Function
    /// \param x x value
    /// \return sigm(x)
    inline double static sigm(const double &x) {
        const double width = 5;
        return x / (abs(x) + width / 2);
    }

    /// Build A submatrices
    inline void Build_A() {
#pragma unroll
        for (size_t i = 0; i < Gsize - 1; i++) {
            double Mu_loc = mu * sigm(i + param - Gsize / 2);
            m(i, i) = Mu_loc + 2 * t;
            m(i + Gsize, i + Gsize) = -Mu_loc - 2 * t;
            m(i, i + 1) = -t;
            m(i + 1, i) = -t;
            m(i + Gsize, i + Gsize + 1) = t;
            m(i + Gsize + 1, i + Gsize) = t;
            m(i, i + 1 + Gsize) = -Delta;
            m(i + 1, i + Gsize) = Delta;
            m(i + Gsize, i + 1) = Delta;
            m(i + 1 + Gsize, i) = -Delta;
        }
        m(Gsize - 1, Gsize - 1) = mu * sigm(Gsize - 1 + param - Gsize / 2) + 2 * t;
        m(2 * Gsize - 1, 2 * Gsize - 1) = -mu * sigm(Gsize - 1 + param - Gsize / 2) - 2 * t;
    }
};


