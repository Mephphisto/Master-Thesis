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
#ifdef TIME_EVOLUTION
#include "Time_Evolution.hpp"
#endif
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

    const size_t num_steps = OMP_NUM_THREADS;
    Vec omegas(num_steps);
    for (int k = 0; k < num_steps; k++) {
        double aux = std::pow(2.0, 14*(k-.5*num_steps)/num_steps-1);
        omegas[k] = aux;
    }
    Mat out = Do_TE(omegas);
     try {
        std::fstream csv_file("Rho_Decay.csv",
                              std::fstream::out);
        assert(csv_file.is_open());
         csv_file << "{{w, maj1-m1j2, maj1+m1j2,norm},";
        for (auto k = 0; k < omegas.size(); k++) {
            csv_file << "{" << omegas[k] << "," << out(0,k) << "," << out(1,k) << "," << out(2,k) << "},";
        }
        csv_file << "}" ;
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


