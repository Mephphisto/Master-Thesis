//
//  Hamiltonian_Matrix.hpp
//  TST
//
//  Created by Jakob Teuffel on 10.12.19.
//  Copyright © 2019 Jakob Teuffel. All rights reserved.
//

#pragma once
#define DBL_EPSILON __DBL_EPSILON__

#include <Eigen/Sparse>
#include <iostream>

class Hamiltonian_MatrixSparse {
protected:
    size_t Msize;
    Eigen::SparseMatrix<std::complex<double>> m;
public:
    Hamiltonian_MatrixSparse(){}
    Hamiltonian_MatrixSparse(size_t size, double t, double mu, double delta){}
    virtual double trace_A() = 0;
    virtual Eigen::SparseMatrix<std::complex<double>> get() = 0;
    virtual bool verify_hermitiity()
    {
        bool res = false;
        for (size_t k = 0; k < 2 * Msize; k++)
        {
            for (size_t j = k; j < 2 * Msize; j++)
            {
                auto diff = abs( m.coeffRef(j, k) - std::conj( m.coeffRef(k, j)));
                if ( diff > __DBL_EPSILON__)
                {
                    std::cout << "Hermitiity error m( " << k << " , " << j << " ) = " << m.coeffRef(k, j) << " m.coeffRef( " << j << " , " << k << " ) = " << m.coeffRef(j, k) << " diff = " << diff << std::endl;
                    res = true;
                }
            }
        }
        return res;
    }
    virtual bool verify_AMatrices()
    {
        bool res = false;
        for (size_t k = 0; k < Msize; k++)
        {
            for (size_t j = 0; j < Msize; j++)
            {
                auto diff = abs( m.coeffRef(j, k) + std::conj( m.coeffRef(j + Msize, k + Msize)));
                if (diff > __DBL_EPSILON__)
                {
                    std::cout << "A Matrix error  m( " << j << " , " << k << " ) = " << m.coeffRef(j, k) << " m.coeffRef( " << j + Msize << " , " << k + Msize << " ) = " << m.coeffRef(j + Msize, k + Msize) << " diff = " << diff << std::endl;
                    res = true;
                }
            }
        }
        return res;
    }
    virtual bool verify_BMatrices()
    {
        bool res = false;
        for (size_t k = Msize; k < 2 * Msize; k++)
        {
            for (size_t j = 0; j < Msize; j++)
            {
                auto diff = abs( m.coeffRef(j, k) + std::conj( m.coeffRef(j + Msize, k - Msize)));
                if (diff > __DBL_EPSILON__)
                {
                    std::cout << "B Matrix error m( " << j << " , " << k << " ) = " << m.coeffRef(j, k) << "; m.coeffRef( " << j + Msize << " , " << k - Msize << " ) = " << m.coeffRef(j + Msize, k - Msize) << " diff = " << diff << std::endl;
                    res = true;
                }
            }
        }
        return res;
    }
};
