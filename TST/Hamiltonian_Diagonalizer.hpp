//
// Created by Jakob Teuffel on 18.12.19.
//

#pragma once

#include "Parameters.hpp"
#include <type_traits>
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <vector>
#include <chrono>
#include "Hamiltonian_MatrixSparse.hpp"
#include "Hamiltonian_Matrix.hpp"

template <class T>
class Diagonalize_Hamiltonian{
protected:
    Eigen::MatrixXd All_EigenValues;
    std::vector<Eigen::VectorXcd> All_EigenVectors;
    Eigen::VectorXd t_s;

#ifdef DEBUG_ACTIVE
    std::chrono::time_point<std::chrono::system_clock>  start;
#endif
    void Preppare() {
#ifdef SHOW_MATRIX_AND_QUIT
        //Build Matrix shom Matrix in CMD and Quit
        T M = T(MATRIX_SIZE, T_COUPLE ,T_START, MU, DELTA);
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
            auto Mat = T(MATRIX_SIZE, T_COUPLE , T_START, MU, DELTA).get();
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
    return;
#endif
#ifdef MU_TO_CSV
        try {
            std::fstream csv_file("MU_M" + std::to_string(MATRIX_SIZE)
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
            Eigen::MatrixXcd M = T(MATRIX_SIZE, T_COUPLE, T_START, MU, DELTA).get();
            for (auto a : M.diagonal()) {
                csv_file << "," << a;
            }
            csv_file << std::endl;


            csv_file.close();
        } catch (...) {
            std::cout << "Error writing MU CSV" << std::endl;
        }
#endif
#ifdef DELTA_TO_CSV
        try {
            std::fstream csv_file("DELTA_M" + std::to_string(MATRIX_SIZE)
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
            Eigen::MatrixXcd M = T(MATRIX_SIZE, T_COUPLE, T_START, MU, DELTA).get();
            for (auto a : M.diagonal<MATRIX_SIZE / 2 + 1>()) {
                csv_file << "," << a;
            }
            for (auto a : M.diagonal<MATRIX_SIZE / 2 - 1>()) {
                csv_file << "," << a;
            }
            std::cout << "MATRIX_SIZE = " << MATRIX_SIZE << " Delta size " << M.diagonal<MATRIX_SIZE / 2 - 1>().size() +
                                                                                M.diagonal<MATRIX_SIZE / 2 + 1>().size() << std::endl;
            csv_file << std::endl;


            csv_file.close();
        } catch (...) {
            std::cout << "Error writing MU CSV" << std::endl;
        }
#endif

#ifdef DEBUG_ACTIVE
        // get start time for runtime Calculation
        this-> start = std::chrono::system_clock::now();
#endif

    }

    virtual void Compute() = 0;

    void Resolve() {
#ifdef EVAL_BY_CMD
        //print Eigenvalues to CMD
            std::cout << "eigenvalues:" << std::endl << All_EigenValues << std::endl;
#endif


#ifdef EVAL_EVAL_BY_CSV
        //writes Eigenvalues to CSV - File for
        try {
            std::fstream csv_file("EigenValues_M" + std::to_string(MATRIX_SIZE)
                                  + "_Tres" + std::to_string(T_RES) + ".csv",
                                  std::fstream::out);
            assert(csv_file.is_open());
            //Write Headder
            csv_file << "Eigenvalues of " << std::to_string(MATRIX_SIZE) << "x"
                     << std::to_string(MATRIX_SIZE) << "Matrix, "
                     << "for t = [" << std::to_string(T_START) << "," << std::to_string(T_END) << "] in "
                     << std::to_string(T_RES) << " Steps. With: "
                     << "t=  " << std::to_string(MU) << "; Delta = " << std::to_string(DELTA)
                     << "using " << std::to_string(OMP_NUM_THREADS) << " OpenMP Threads" << std::endl;
            csv_file << "M" << "," << std::to_string(MATRIX_SIZE) << " Eigenvalues  ... " << std::endl;
            //write Eigenvalues
            for (auto k = 0; k < T_RES; k++) {
                csv_file << t_s(k);
                for (auto a : All_EigenValues.col(k)) {
                    csv_file << "," << a;
                }
                csv_file << std::endl;
            }
            csv_file.close();
        } catch (...) {
            std::cout << "Error writing EigenValues CSV" << std::endl;
        }
#endif
#ifdef EVAL_EVEC_BY_CSV
        try {
            std::fstream csv_file("EigenVectors_M" + std::to_string(MATRIX_SIZE)
                                  + "_Tres" + std::to_string(T_RES) + ".csv",
                                  std::fstream::out);
            assert(csv_file.is_open());
            //Write Headder
            csv_file << "Eigenvalues of " << std::to_string(MATRIX_SIZE) << "x" << std::to_string(MATRIX_SIZE)
                     << "Matrix, "
                     << "for t = [" << std::to_string(T_START) << "," << std::to_string(T_END) << "] in "
                     << std::to_string(T_RES) << " Steps. With: "
                     << "t=  " << std::to_string(MU) << "; Delta = " << std::to_string(DELTA)
                     << "using " << std::to_string(OMP_NUM_THREADS) << " OpenMP Threads" << std::endl;
            csv_file << "M" << "," << std::to_string(MATRIX_SIZE) << " Eigenvalues  ... " << std::endl;
            //write Eigenvalues
            std::cout << "Min Eval " << All_EigenValues.minCoeff() << std::endl;
            for (auto V : All_EigenVectors) {
                for (auto a : V) {
                    csv_file << "," << a;
                }
                csv_file << std::endl;
            }

            csv_file.close();
        } catch (...) {
            std::cout << "Error writing EigenVectors CSV" << std::endl;
        }

#endif

#ifdef DEBUG_ACTIVE
        //Compute Calculation Time adn print to CMD
        std::cout << std::endl << std::endl << "Runtime = "
                  <<
                  std::chrono::duration_cast<std::chrono::milliseconds>(
                          std::chrono::system_clock::now() - start).count()
                  << "ms" <<
                  std::endl;
#endif

#if defined EVAL_BY_CMD || defined DEBUG_ACTIVE
        //Wait for ENTER before Closing window
        std::cout << "FINISHED!" << std::endl;
#endif
    }

public:
    Diagonalize_Hamiltonian() {
        bool right_T = std::is_base_of<Hamiltonian_Matrix, T>::value;
        right_T |= std::is_base_of<Hamiltonian_MatrixSparse, T>::value;
        assert(right_T);
    }
     void Do() {
        All_EigenValues.resize(MATRIX_SIZE,T_RES);
        t_s.resize(T_RES);
        Preppare();
        this->Compute();
        Resolve();
    }
    virtual ~Diagonalize_Hamiltonian() {}
};