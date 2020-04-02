


#pragma once

#include "Hamiltonian_Matrix.hpp"
#include <Eigen/Dense>

class _1DTopSuperConMatrix : Hamiltonian_Matrix {
public:
    _1DTopSuperConMatrix(size_t size_in, double t_con, double param, double mu_in, double Delta_in) {
        this->t = -t_con;
        this->mu = -mu_in;
        this->Delta = Delta_in / 2;
        this->param = param;
        Gsize = size_in / 2;
        m = Eigen::MatrixXcd(size_in, size_in);
        m.setZero();
        Build_A();
        //Build_B();
    }

    Eigen::MatrixXcd get() override {
        return m;
    }

    double trace_A() override {
        return 0;
    }

private:
    Eigen::MatrixXcd m;
    double t, mu, Delta, param;
    size_t Gsize;
    inline double sigm(const  double & x){
        const double width = 5;
        return x /(abs(x) + width / 2);
    }
    void Build_A() {
        for (size_t i = 0; i < Gsize - 1; i++) {
            double Mu_loc = mu * sigm(i + param - Gsize / 2);
            m(i, i) = Mu_loc + 2 * t;
            m(i + Gsize, i + Gsize) = -Mu_loc - 2*t;
            m(i, i + 1) = -t;
            m(i + 1, i) = -t;
            m(i + Gsize, i + Gsize + 1) = t;
            m(i + Gsize + 1, i + Gsize) = t;
            m(i, i + 1 + Gsize) = -Delta;
            m(i + 1, i + Gsize) = Delta;
            m(i + Gsize, i + 1) = Delta;
            m(i + 1 + Gsize, i) = -Delta;
        }
        m(Gsize - 1, Gsize - 1) = mu*sigm(Gsize-1+param-Gsize/2) + 2*t;
        m(2 * Gsize - 1, 2 * Gsize - 1) = -mu*sigm(Gsize-1+param-Gsize/2) - 2*t;
    }
};


