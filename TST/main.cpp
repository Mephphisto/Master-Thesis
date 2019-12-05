//
//  main.cpp
//  TST
//
//  Created by Jakob Teuffel on 05.12.19.
//  Copyright © 2019 Jakob Teuffel. All rights reserved.
//
#define EIGEN_USE_MKL

#include <iostream>
#include <Eigen/Dense>
using Eigen::MatrixXd;
int main()
{
    MatrixXd m(2,2);
    m(0,0) = 0;
    m(1,0) = 1;
    m(0,1) = -1;
    m(1,1) = 0;
    std::cout << m << std::endl;
}
