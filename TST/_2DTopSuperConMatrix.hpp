//
// Created by Jakob Teuffel on 16.12.19.
//


#pragma once
#include "Hamiltonian_Matrix.hpp"
#include <complex>
#include <Eigen/Dense>

class _2DTopSuperConMatrix :  Hamiltonian_Matrix
{
private:
    Eigen::MatrixXcd m;
    double t, mu, Delta;
    size_t Gsize , Msize;

    void Build_A()
    {
        for(size_t k = 0 ; k <  Msize ; k++){
            for (size_t j = 0 ; j <  Msize ; j++){
                size_t x_j = index_x(j);
                size_t y_j = index_y(j);
                size_t x_k = index_x(k);
                size_t y_k = index_y(k);
                if ((x_j == x_k) && (y_j == y_k)) {
                    m(k, j) = -mu;
                    m(k + Msize, j + Msize) = -mu;
                } else if (( x_k == x_j) && (y_k == x_k +1)) {
                    m(k, j) = -t;
                    m(k + Msize, j + Msize) = -t;
                } else if (( x_k == x_j + 1) && (y_k == x_k)){
                    m(k, j) = -t;
                    m(k + Msize, j + Msize) = -t;
                } else {
                    m( k, j)= 0;
                    m( k + Msize, j + Msize);
                }
            }
        }
    }
    void Build_B(){
        for( size_t k = 0 ; k <  Msize ; k++){
            for ( size_t j = Msize ; j < 2 * Msize ; j++){
                size_t x_j = index_x(j);
                size_t y_j = index_y(j);
                size_t x_k = index_x(k);
                size_t y_k = index_y(k);
                if  (( x_k == x_j) && (y_k == x_k +1)){
                    m(k, j) = Delta;
                    m(k + Msize, j - Msize) = std::conj(Delta);
                } else if (( x_k == x_j + 1) && (y_k == x_k)){
                    m(k, j) = - std::conj(Delta);
                    m(k + Msize, j - Msize) = -Delta;
                } else {
                    m(k,j)= 0;
                }
            }
        }
    }
    inline size_t at( size_t index_x , size_t index_y){
        return  index_x + Gsize * index_y;
    }
    inline size_t index_x (size_t pos){
        return pos % Gsize;
    }
    inline size_t index_y (size_t pos){
        return pos / Gsize;
    }
public:
    _2DTopSuperConMatrix(size_t size_in, double t_in, double mu_in, double Delta_in)
    {
        this -> t = -t_in;
        this -> mu = -mu_in;
        this -> Delta = Delta_in / 2;
        this -> Msize = size_in / 2;
        this -> Gsize = sqrt(Msize);

        m = Eigen::MatrixXcd(size_in, size_in);
        Build_A();
        Build_B();
    }
    Eigen::MatrixXcd get() override
    {
        return m;
    }
    double trace_A() override
    {
        return 0;
    }
};

