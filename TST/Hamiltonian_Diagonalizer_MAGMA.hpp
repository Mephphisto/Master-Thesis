//
// Created by Jakob Teuffel on 18.12.19.
//

#pragma once

#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include "TST_Magma_Solver.hpp"
#include <magma_lapack.h>
#include <chrono>

template<class T>
void Diagonalize_Hamiltonian_magma() {
#ifdef SHOW_MATRIX_AND_QUIT
    //Build Matrix shom Matrix in CMD and Quit
    T M = T(MATRIX_SIZE, T_START, MU, DELTA);
    bool Mat_Error = false;
    Mat_Error = M.verify_hermitiity() || Mat_Error;
    Mat_Error = M.verify_AMatrices() || Mat_Error;
    Mat_Error|= M.verify_BMatrices() || Mat_Error;
    if (!Mat_Error)
    {
        std::cout << "No Matrix Errors found" << std::endl;
    }else
    {
        std::cin.ignore();
    }
    std::cout << M.get();
    std::cin.get();
    return;
#endif
#ifdef MATRIX_TO_CSV
    {
        try {
            std::fstream csv_file("Matrix" + std::to_string(MATRIX_SIZE)
                                  + "_Tres" + std::to_string(T_RES) + ".csv",
                                  std::fstream::out);
            assert(csv_file.is_open());
            //Write Headder
            csv_file << "Matrix of " << std::to_string(MATRIX_SIZE) << "x" << std::to_string(MATRIX_SIZE) << "Matrix, "
                     << "for t = " << std::to_string(T_START) << " With: "
                     << "MU=  " << std::to_string(MU) << "; Delta = " << std::to_string(DELTA)
                     << std::endl;
            csv_file << "M" << "," << std::to_string(MATRIX_SIZE) << " Eigenvalues  ... " << std::endl;
            //write Matrix
            auto Mat = T(MATRIX_SIZE, T_START, MU, DELTA).get();
            for (auto k = 0; k < MATRIX_SIZE; k++) {
                for (auto a : Mat.col(k)) {
                    csv_file << "," << a;
                }
                csv_file << std::endl;
            }
            csv_file.close();
        } catch (...) {
            std::cout << "Error Closing file" << std::endl;
        }
    }
#endif

#ifdef DEBUG_ACTIVE
    // get start time for runtime Calculation
    auto start = std::chrono::system_clock::now();

    // print number of used threads
    std::cout << "Running " << Eigen::nbThreads() << " Eigen Threads" << std::endl;
#endif

    // Storage for Computed Eigen Values
    Eigen::MatrixXd All_EigenValues(MATRIX_SIZE, T_RES);
    std::vector<Eigen::MatrixXcd> All_EigenVectors(T_RES);
    Eigen::VectorXd t_s(T_RES);
    {
        // Create Solver
        magma::SelfAdjointEigenSolver<Eigen::MatrixXcd> Solver(MATRIX_SIZE);
        for (auto k = 0; k < T_RES; k++) {
#ifdef DEBUG_ACTIVE
            //verify accurate stride through sample space
            {
                //Build String first and print then Prevents Race Condition!
                std::string msg;
                msg = " Status= " + std::to_string(static_cast<float>(k) / T_RES * 100) + "%\n";
                std::cout << msg;
            }
#endif
            //Storage for Matrix and Trace
            Eigen::MatrixXcd m;
            double tr;
            double t = T_START + (T_END - T_START) * (double(k) / double(T_RES - 1));
            {
                //Build Matrix to be Solved by MKL
                T M = T(MATRIX_SIZE, t, MU, DELTA);
                m = M.get();
                tr = M.trace_A();
            }
            //Solve the Matrix
            Solver.compute(m);
            // Fetch  Eigenvalues from Solver
            Eigen::VectorXd EigenValues = Solver.eigenvalues().col(0).real();

            //Save Eigenvalues and add constant Correction terms
            All_EigenValues.col(k) = EigenValues + Eigen::VectorXd::Constant(MATRIX_SIZE, (tr - EigenValues.sum()) / 2);
            All_EigenVectors[k] = EigenVectors;
            t_s(k) = t;
        }
    }


#ifdef DEBUG_ACTIVE
    // Save End Time
    auto end = std::chrono::system_clock::now();
#endif

#ifdef EVAL_BY_CMD
    //print Eigenvalues to CMD
        std::cout << "eigenvalues:" << std::endl << All_EigenValues << std::endl;
#endif


#ifdef    EVAL_BY_CSV
    //writes Eigenvalues to CSV - File for
    try {
        std::fstream csv_file("EigenValues_M" + std::to_string(MATRIX_SIZE)
                              + "_Tres" + std::to_string(T_RES) + ".csv",
                              std::fstream::out);
        assert(csv_file.is_open());
        //Write Headder
        csv_file << "Eigenvalues of " << std::to_string(MATRIX_SIZE) << "x" << std::to_string(MATRIX_SIZE) << "Matrix, "
                 << "for t = [" << std::to_string(T_START) << "," << std::to_string(T_END) << "] in "
                 << std::to_string(T_RES) << " Steps. With: "
                 << "t=  " << std::to_string(MU) << "; Delta = " << std::to_string(DELTA)
                 << "using " << std::to_string(OMP_NUM_THREADS) << " OpenMP Threads" << std::endl;
        csv_file << "M" << "," << std::to_string(MATRIX_SIZE) << " Eigenvalues  ... " << std::endl;
        //write Eigenvalues
        for (auto k = 0; k <= T_RES; k++) {
            csv_file << t_s(k);
            for (auto a : All_EigenValues.col(k)) {
                csv_file << "," << a;
            }
            csv_file << std::endl;
        }
        csv_file.close();
    }catch (...){
        std::cout << "Error writing EigenValues CSV"  << std::endl;
    }
    try {
        std::fstream csv_file("EigenVectors_M" + std::to_string(MATRIX_SIZE)
                              + "_Tres" + std::to_string(T_RES) + ".csv",
                              std::fstream::out);
        assert(csv_file.is_open());
        //Write Headder
        csv_file << "Eigenvalues of " << std::to_string(MATRIX_SIZE) << "x" << std::to_string(MATRIX_SIZE) << "Matrix, "
                 << "for t = [" << std::to_string(T_START) << "," << std::to_string(T_END) << "] in "
                 << std::to_string(T_RES) << " Steps. With: "
                 << "t=  " << std::to_string(MU) << "; Delta = " << std::to_string(DELTA)
                 << "using " << std::to_string(OMP_NUM_THREADS) << " OpenMP Threads" << std::endl;
        csv_file << "M" << "," << std::to_string(MATRIX_SIZE) << " Eigenvalues  ... " << std::endl;
        //write Eigenvalues
        for (auto k = 0; k <= T_RES; k++) {
            for (auto a : All_EigenVectors) {
                csv_file << "," << a;
            }
            csv_file << std::endl;
        }
        csv_file.close();
    }catch (...){
        std::cout << "Error writing EigenVectors CSV"  << std::endl;
    }
#endif

#ifdef DEBUG_ACTIVE
    //Compute Calculation Time adn print to CMD
    std::cout << std::endl << std::endl << "Runtime = "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << "ms" << std::endl;
#endif

#if defined EVAL_BY_CMD || defined DEBUG_ACTIVE
    //Wait for ENTER before Closing window
    std::cout << "FINISHED!" << std::endl;
#endif
}

