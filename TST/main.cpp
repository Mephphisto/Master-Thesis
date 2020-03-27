//
//  main.cpp
//  TST
//
//  Created by Jakob Teuffel on 05.12.19.
//  Copyright © 2019 Jakob Teuffel. All rights reserved.
//  Copyright © 2019 Jakob Teuffel. All rights reserved.
//


#include "Parameters.hpp"
#include "_1DTopSuperConMatrix.hpp"
#include "_2DTopSuperConMatrix.hpp"
#include "_2DTopSuperConMatrixSparse.hpp"


#ifdef USE_MAGMA
#include "Hamiltonian_Diagonalizer_MAGMA.hpp"
#elif defined USE_OpenCL
#include "Hamiltonian_Diagonalizer_OpenCl.hpp"
#else
#include "Hamiltonian_Diagonalizer_Eigen.hpp"
#endif

int main() {
#ifdef USE_MAGMA
    Diagonalize_Hamiltonian_magma<_1DTopSuperConMatrix>().Do();
#elif defined USE_OpenCL
    Diagonalize_Hamiltonian_OpenCL<_2DTopSuperConMatrix>();
#else
    std::cout << "Num Eingen Threads " << Eigen::nbThreads( ) << std::endl;
    Diagonalize_Hamiltonian_Eigen<_1DTopSuperConMatrix>().Do();
#endif
    return 0;
}