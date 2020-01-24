//
//  Hamiltonian_Matrix.hpp
//  TST
//
//  Created by Jakob Teuffel on 10.12.19.
//  Copyright © 2019 Jakob Teuffel. All rights reserved.
//

#ifndef Hamiltonian_Matrix_hpp
#define Hamiltonian_Matrix_hpp
#define DBL_EPSILON __DBL_EPSILON__

#include <Eigen/Dense>
#include <iostream>

class Hamiltonian_Matrix {
protected:
    size_t Msize;
    Eigen::MatrixXcd m;
public:
    Hamiltonian_Matrix(){}
    Hamiltonian_Matrix(size_t size, double t, double mu, double delta){};
    virtual double trace_A() = 0;
    virtual Eigen::MatrixXcd get() = 0;
    virtual bool verify_hermitiity()
    {
        bool res = false;
        for (size_t k = 0; k < 2 * Msize; k++)
        {
            for (size_t j = k; j < 2 * Msize; j++)
            {
                auto diff = abs(m(j, k) - conj(m(k, j)));
                if ( diff > DBL_EPSILON)
                {
                    std::cout << "Hermitiity error m( " << k << " , " << j << " ) = " << m(k, j) << " m( " << j << " , " << k << " ) = " << m(j, k) << " diff = " << diff << std::endl;
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
                auto diff = abs(m(j, k) + conj(m(j + Msize, k + Msize)));
                if (diff > DBL_EPSILON)
                {
                    std::cout << "A Matrix error  m( " << j << " , " << k << " ) = " << m(j, k) << " m( " << j + Msize << " , " << k + Msize << " ) = " << m(j + Msize, k + Msize) << " diff = " << diff << std::endl;
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
                auto diff = abs(m(j, k) + conj(m(j + Msize, k - Msize)));
                if (diff > DBL_EPSILON)
                {
                    std::cout << "B Matrix error m( " << j << " , " << k << " ) = " << m(j, k) << "; m( " << j + Msize << " , " << k - Msize << " ) = " << m(j + Msize, k - Msize) << " diff = " << diff << std::endl;
                    res = true;
                }
            }
        }
        return res;
    }
};

#endif /* Hamiltonian_Matrix_h */
