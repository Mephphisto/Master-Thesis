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

class Hamiltonian_Matrix {
public:
    Hamiltonian_Matrix(){}
    Hamiltonian_Matrix(size_t size, double t, double mu, double delta){};
    virtual double trace_A() = 0;
    virtual Eigen::MatrixXcd get() = 0;
};

#endif /* Hamiltonian_Matrix_h */
