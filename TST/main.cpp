//
//  main.cpp
//  TST
//
//  Created by Jakob Teuffel on 05.12.19.
//  Copyright © 2019 Jakob Teuffel. All rights reserved.
//
#define EIGEN_USE_MKL
#define MATRIX_SIZE 1000

#include <iostream>
#include <Eigen/Dense>
#include <omp.h>

int main()
{
    omp_set_num_threads(4);
    Eigen::MatrixX<std::complex<double>>m(MATRIX_SIZE,MATRIX_SIZE);
    #pragma omp parallel for
    for (size_t i = 0; i < MATRIX_SIZE; i++) {
        #pragma omp parallel for
        for (size_t j = 0; j < MATRIX_SIZE; j++) {
            double z = double(j)-double(i);
            double x = double(j)+double(i);
            m(i,j) = std::complex<double>(1+(x-MATRIX_SIZE/2-2)/(z*z+1),0)*exp(std::complex<double>(0,z));
        }
    }
    //std::cout << m << std::endl;
    Eigen::ComplexEigenSolver<Eigen::MatrixX<std::complex<double>>> Solver(MATRIX_SIZE);
    Solver.compute(m);
    Eigen::VectorX<double> EigenValues = Solver.eigenvalues().col(0).real();
    std::sort(EigenValues.data(),EigenValues.data()+EigenValues.size());
    std::cout << "eigenvalues:" << std::endl << EigenValues << std::endl;
}
