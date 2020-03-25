


#pragma once

#include "Hamiltonian_Matrix.hpp"
#include <Eigen/Dense>

class _1DTopSuperConMatrix : public Hamiltonian_Matrix {
public:
    _1DTopSuperConMatrix(size_t size_in, double t_coup, double t_in, double mu_in, double Delta_in) {
        this->t = -t_in;
        this->t_coup = t_coup;
        this->mu = -mu_in;
        this->Delta = Delta_in / 2;
        Msize = size_in / 2;
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
    double t, mu, Delta, t_coup;

    void Build_A() {
        for (size_t i = 0; i < Msize - 1; i++) {
            m(i, i) = -mu;
            if (i = (size_t) Msize / 2) {
                m(i, i + 1) = -t;
                m(i + 1, i) = -t;
                m(i + Msize, i + Msize + 1) = t;
                m(i + Msize + 1, i + Msize) = t;
            } else if (i < Msize / 2) {
                m(i, i + 1) = -t_coup;
                m(i + 1, i) = -t_coup;
                m(i + Msize, i + Msize + 1) = t_coup;
                m(i + Msize + 1, i + Msize) = t_coup;
            } else {
                m(i, i + 1) = t_coup;
                m(i + 1, i) = t_coup;
                m(i + Msize, i + Msize + 1) = -t_coup;
                m(i + Msize + 1, i + Msize) = -t_coup;
            }
            m(i + Msize, i + Msize) = mu;
            m(i, i + 1 + Msize) = -Delta;
            m(i + 1, i + Msize) = Delta;
            m(i + Msize, i + 1) = Delta;
            m(i + 1 + Msize, i) = -Delta;
        }
        m(Msize - 1, Msize - 1) = -mu;
        m(2 * Msize - 1, 2 * Msize - 1) = mu;
    }
};

