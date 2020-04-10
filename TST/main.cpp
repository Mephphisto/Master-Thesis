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
    Vec omegas (48*4);
    for (size_t k = -96; k < 96; k++){
        omegas[k+96] = pow(2, k);
    }
    Vec out = Do(omegas);
    try {
            std::fstream csv_file("Rho_Decay.csv",
            std::fstream::out);
            assert(csv_file.is_open());
            for (auto k = 0; k < 48*4; k++) {
                csv_file << "{" << omegas[k] << "," << out[k] << "},";
                csv_file << std::endl;
            }
            csv_file.close();
    } catch (...) {
        std::cout << "Error writing EigenValues CSV" << std::endl;
    }
#endif
    return 0;
}


