//
//  Hamiltonian_Matrix.hpp
//  TST
//
//  Created by Jakob Teuffel on 10.12.19.
//  Copyright © 2019 Jakob Teuffel. All rights reserved.
//

#ifndef Hamiltonian_Matrix_hpp
#define Hamiltonian_Matrix_hpp

#include <Eigen/Dense>
#include <iostream>
#include "Typedefs.hpp"

/// Abstract Base class that holds the Hermitian Matrix Mat_cd representig our Hamiltonian.
/// The Matrix has the form
/// \f[
/// \left(
/// \begin{array}{cc}
/// A & B \\
/// \text{$\epsilon $B}^* & \text{$\epsilon $A}^* \\
/// \end{array}
/// \right) ;
///    \epsilon =
///  \begin{cases}
/// 1 & boson \\
/// -1 & fermion
/// \end{cases}
/// \f]
/// coontains Interface and functions to check Hermiticity.
class Hamiltonian_Matrix {
protected:
    /// Size - n of the n x n Matrix
    size_t Msize;
    /// The "Eigen" Matrix that holds the entries.
    Mat_cd m;
public:
    /// Constructor
    /// \return Hamiltonian_Matrix Opbject need to extract Eigen Mat_cd with "get()"
    Hamiltonian_Matrix(const size_t &size) {
        this->Msize = size / 2;
        m = Mat_cd::Zero(size, size);
    };

    /// All derived classes need a trace funtion as the trace represents a necceccary
    /// Energy shift. from the commutator Relation \f$c^+_j  c_i = 1/2 ( c^+_i c_i + c_i c^+_i - [c_i, c^+_i])\f$
    virtual double trace_A() = 0;

    /// All derived classes need to be able to return the Eigen Storage Mat_cd for computations
    virtual Mat_cd &get() = 0;

    /// This Function checks Hermiticity fot all derived Classes
    /// \return True if Matrix is Non Hermitian
    virtual bool verify_hermitiity() {
        bool res = false; /// Result of chek false coresponds to a corect hermitian Matrix
        for (size_t k = 0; k < 2 * Msize; k++) {
            for (size_t j = k; j < 2 * Msize; j++) {
                auto diff = std::abs(m(j, k) - std::conj(m(k, j)));
                if (diff > __DBL_EPSILON__) {
                    std::cout << "Hermitiity error m( " << k << " , " << j << " ) = " << m(k, j) << " m( " << j << " , "
                              << k << " ) = " << m(j, k) << " diff = " << diff << std::endl;
                    res = true;
                }
            }
        }
        return res;
    }

    /// Verifies, that the the \f$ A \f$ and \f$ \epsilon A^* \f$ submatrices are indeed related as desired.
    /// \return True if Mat_cd is incorrect
    virtual bool verify_AMatrices() {

        bool res = false;
        for (size_t k = 0; k < Msize; k++) {
            for (size_t j = 0; j < Msize; j++) {
                auto diff = std::abs(m(j, k) + std::conj(m(j + Msize, k + Msize)));
                if (diff > __DBL_EPSILON__) {
                    std::cout << "A Matrix error  m( " << j << " , " << k << " ) = " << m(j, k) << " m( " << j + Msize
                              << " , " << k + Msize << " ) = " << m(j + Msize, k + Msize) << " diff = " << diff
                              << std::endl;
                    res = true;
                }
            }
        }
        return res;
    }

    /// Verifies, that the the \f$ B \f$ and \f$ \epsilon B^* \f$ submatrices are indeed related as desired.
    /// \return True if Mat_cd is incorrect
    virtual bool verify_BMatrices() {
        bool res = false;
        for (size_t k = Msize; k < 2 * Msize; k++) {
            for (size_t j = 0; j < Msize; j++) {
                auto diff = std::abs(m(j, k) + std::conj(m(j + Msize, k - Msize)));
                if (diff > __DBL_EPSILON__) {
                    std::cout << "B Matrix error m( " << j << " , " << k << " ) = " << m(j, k) << "; m( " << j + Msize
                              << " , " << k - Msize << " ) = " << m(j + Msize, k - Msize) << " diff = " << diff
                              << std::endl;
                    res = true;
                }
            }
        }
        return res;
    }
};

#endif /* Hamiltonian_Matrix_h */
