//
// Created by jakob on 15.10.20.
//

#ifndef TST_MAJORANA_TRANSFORM_HPP
#define TST_MAJORANA_TRANSFORM_HPP

#include <Eigen/Dense>
#include <iostream>
#include <omp.h>
#include "Typedefs.hpp"

inline cd MyDot(Vec_cd X, Vec_cd Y) {
    cd res1 = 0, res2 = 0;
    size_t L = X.size() / 2;
    assert(Y.size() == 2 * L);
    for (size_t k = 0; k < L / 2; k++) {
        cd norm1 = X[k] * Y[k];
        res1 += norm1;
    }
    for (size_t k = L; k < 3 * L / 2; k++) {
        cd norm2 = X[k] * Y[k];
        res1 += norm2;
    }
    cd res = res1 + res2;
    return res;
}

inline std::tuple<Vec_cd, Vec_cd> Majoranaize(Vec_cd v, Vec_cd w) {
    Vec_cd X, Y;
    {

        {
            X = v - (MyDot(v, w) / MyDot(w, w)) * w;
            X.normalize();
        }
        {
            cd dotXX = X.transpose() * X;
            cd dotvX = (v.transpose() * X);
            cd dotwX = (w.transpose() * X);
            Y = v - dotvX / dotXX * X + w - dotwX / dotXX * X;
            Y.normalize();
        }
#ifdef DEBUG_ACTIVE
#pragma omp critical
        {
            double norm1 = std::abs(MyDot(X.conjugate(), X)), norm2 = std::abs(MyDot(Y.conjugate(), Y));
            if (static_cast<int>(std::log(norm1)) >= static_cast<int>(std::log(norm1))) {
                std::cerr <<  " For the Majorana Transform to be successfull |maj1_s|² << |maj2_s|² is neede ";
            }

            std::cout << "|maj1_s|²= " << norm1
                      << " |maj2_s|²= " << norm2
                      << std::endl;
        }
    }
#endif
    return std::make_tuple(X, Y);
}

}
#endif //TST_MAJORANA_TRANSFORM_HPP
