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
#include "Time_Evolution.hpp"

#ifdef USE_MAGMA

#include "Hamiltonian_Diagonalizer_MAGMA.hpp"

#elif defined USE_OpenCL
#include "Hamiltonian_Diagonalizer_OpenCl.hpp"
#else
#include "Hamiltonian_Diagonalizer_Eigen.hpp"
#include "Hamiltonian_Diagonalizer_PLASMA.hpp"
#endif

int main() {
#ifndef TIME_EVOLUTION
#ifdef USE_MAGMA
    Diagonalize_Hamiltonian_magma<HAMILTONIAN>().Do();
#elif defined USE_OpenCL
    Diagonalize_Hamiltonian_OpenCL<HAMILTONIAN>().Do();
#else
    Diagonalize_Hamiltonian_Eigen<HAMILTONIAN>().Do();
#endif
#else
    Vec omegas(5);
    omegas << 0.01, 1.0, 2.0, 3.0, 4.0;
    Vec out = Do(omegas);
    std::cout << out;
#endif
    return 0;
}


