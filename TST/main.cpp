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
#include <math.h>


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
    std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
    const size_t num_steps = 48;
    Vec omegas(num_steps * 4);
    for (int k = 0; k < num_steps * 4; k++) {
        double aux = std::pow(2.0, k-2.0*num_steps);
        omegas[k] = aux;
    }
    Vec out = Do(omegas);
    try {
        std::fstream csv_file("Rho_Decay.csv",
                              std::fstream::out);
        assert(csv_file.is_open());
        for (auto k = 0; k < 4 * num_steps; k++) {
            csv_file << "{" << omegas[k] << "," << out[k] << "},";
            csv_file << std::endl;
        }
        csv_file.close();
    } catch (...) {
        std::cout << "Error writing EigenValues CSV" << std::endl;
    }
    std::cout << std::endl << std::endl << "Runtime = "
              <<
              std::chrono::duration_cast<std::chrono::milliseconds>(
                      std::chrono::system_clock::now() - start).count()
              << "ms" <<
              std::endl;
#endif
    return 0;
}


