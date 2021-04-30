#pragma once

#define OPTIM_ENABLE_EIGEN_WRAPPERS

#include "optim.hpp"
#include <autodiff/forward.hpp>
#include <autodiff/forward/eigen.hpp>

namespace mjize {
    typedef std::complex<double> cd;
    typedef std::complex<autodiff::dual> cdual;

    Eigen::Matrix2<cdual> R, Id;
    Eigen::Matrix4<cdual> R4, Id4;
    cd i = cd(0, 1);
    cdual id = cdual(0, 1);

    autodiff::dual opt_fnd(const autodiff::VectorXdual &x) {
        Eigen::Matrix2<cdual> A = Eigen::Matrix2<cdual>::Zero();
        A << id * x[0], x[1] + id * x[2], cdual(-x[1]) + id * x[2], id * x[3];
        A = autodiff::forward::exp(cdual(0, x[4])) * ((A - Id) * (Id + A).inverse()).eval();
        cdual res1 = cdual(-1) * (A.transpose() * R * A).real().trace();
        autodiff::dual res;
        res.val = res1.real().val;
        res.grad = res1.real().grad;

        return res;
    }

    double opt_fn(const Eigen::VectorXd &x, Eigen::VectorXd *grad_out, void *opt_data) {
        autodiff::dual u;
        autodiff::VectorXdual xd = x.eval();

        if (grad_out) {

            *grad_out = autodiff::forward::gradient(opt_fnd, autodiff::wrt(xd), autodiff::forward::at(xd), u);
        } else {
            u = opt_fnd(xd);
        }

        return u.val;
    }

    autodiff::dual opt_fnd4(const autodiff::VectorXdual &x) {
        Eigen::Matrix4<cdual> A = Eigen::Matrix4<cdual>::Zero();
        //@formatter:off
        A << cdual(0,     x[0]), cdual(x[1],   x[2]),  cdual(x[3],   x[4]),  cdual(x[5],  x[6]),
             cdual(-x[1], x[2]), cdual(0,      x[7]),  cdual(x[8],   x[9]),  cdual(x[10], x[11]),
             cdual(-x[3], x[4]), cdual(-x[8],  x[9]),  cdual(0,      x[12]), cdual(x[13], x[14]),
             cdual(-x[5], x[6]), cdual(-x[10], x[11]), cdual(-x[13], x[14]), cdual(0,     x[15]);
        //@formatter:on
        A = autodiff::forward::exp(cdual(0, x[16])) * ((A - Id4) * (Id4 + A).inverse()).eval();
        cdual res1 = cdual(-1) * (A.transpose() * R4 * A).real().trace();
        autodiff::dual res;
        res.val = res1.real().val;
        res.grad = res1.real().grad;

        return res;
    }

    double opt_fn4(const Eigen::VectorXd &x, Eigen::VectorXd *grad_out, void *opt_data) {
        autodiff::dual u;
        autodiff::VectorXdual xd = x.eval();

        if (grad_out) {

            *grad_out = autodiff::forward::gradient(opt_fnd4, autodiff::wrt(xd), autodiff::forward::at(xd), u);
        } else {
            u = opt_fnd(xd);
        }

        return u.val;
    }

/// Swaps the u and v subvecors thus the first with the last half of a vector.
/// equivalent to gamma = (sigma^x tensor 1_N)
/// \param X Vector
/// \return gamma * Vector
    Vec_cd GammaSwap(const Vec_cd &X) {
        const size_t N = X.size() / 2;
        Vec_cd gX(2 * N);
        gX.head(N) = X.tail(N);
        gX.tail(N) = X.head(N);
        return gX;
    }

    Eigen::MatrixXcd de_dual(const Eigen::MatrixX<cdual> &x) {
        Eigen::MatrixXcd y(x.rows(), x.cols());
        for (int k = 0; k < x.rows(); ++k) {
            for (int j = 0; j < x.cols(); ++j) {
                y(k, j) = cd(x(k, j).real().val, x(k, j).imag().val);
            }
        }
        return y;
    }

    std::tuple<Vec_cd, Vec_cd> Majoranize(const Vec_cd &X, const Vec_cd &Y) {
        Id << cdual(1), cdual(0), cdual(0), cdual(1);
        Vec_cd Xg = GammaSwap(X), Yg = GammaSwap(Y);
        R << cdual(X.cwiseProduct(Xg).sum()), cdual(Y.cwiseProduct(Xg).sum()), cdual(X.cwiseProduct(Yg).sum()), cdual(
                Y.cwiseProduct(Yg).sum());

        Eigen::Matrix2cd r = de_dual(R);
        Eigen::VectorXd params = Eigen::VectorXd::Ones(5);


        optim::algo_settings_t opt_set;
        opt_set.print_level = 3;
        bool success = optim::bfgs(params, opt_fn, nullptr, opt_set);
#ifdef DEBUG_ACTIVE
        if (success) {
            std::cout << "bfgs: autodiff test completed successfully." << std::endl;
        } else {
            std::cout << "bfgs: autodiff test completed unsuccessfully." << std::endl;
        }
#endif
        Eigen::Matrix2cd A = Eigen::Matrix2cd::Zero();
        A << i * params[0], params[1] + i * params[2], -params[1] + i * params[2], i * params[3];
        A = std::exp(cd(0, params[4])) *
                ((A - Eigen::Matrix2cd::Identity()) * (Eigen::Matrix2cd::Identity() + A).inverse()).eval();
#ifdef DEBUG_ACTIVE
        std::cout << "Err_in = " << 4.0 - (r).real().trace() << std::endl;
        std::cout << "Err_out = " << 4.0 - (A.transpose() * r * A).real().trace() << std::endl;
        std::cout << "R  = " << r << std::endl;
        std::cout << "A^T * R * A = " << -A.transpose() * r * A << std::endl;
        std::cout << "A  = " << A << std::endl;
        std::cout << "|X^+ - gamma X|²  = " << (X.conjugate() - GammaSwap(X)).squaredNorm() << std::endl;
        std::cout << "|Y^+ - gamma Y|²  = " << (Y.conjugate() - GammaSwap(Y)).squaredNorm() << std::endl;
        {
            Vec_cd Xpr = (A(0, 0) * X + A(0, 1) * Y),
                    Ypr = (A(1, 0) * X + A(1, 1) * Y);
            std::cout << "|X'^+ - gamma X'|²  = " << (Xpr.conjugate() - GammaSwap(Xpr)).squaredNorm() << std::endl;
            std::cout << "|Y'^+ - gamma Y'|²  = " << (Ypr.conjugate() - GammaSwap(Ypr)).squaredNorm() << std::endl;
        }
#endif
        return std::make_tuple(A(0, 0) * X + A(0, 1) * Y, A(1, 0) * X + A(1, 1) * Y);
    }

    std::tuple<Vec_cd, Vec_cd, Vec_cd, Vec_cd>
    Majoranize(const Vec_cd &X, const Vec_cd &Y, const Vec_cd &Z, const Vec_cd &Q) {
        //@formatter:off
        Id4 <<
            cdual(1), cdual(0), cdual(0), cdual(0),
            cdual(0), cdual(1), cdual(0), cdual(0),
            cdual(0), cdual(0), cdual(1), cdual(0),
            cdual(0), cdual(0), cdual(0), cdual(1);
        //@formatter:on
        {
            Vec_cd Xg = GammaSwap(X), Yg = GammaSwap(Y), Zg = GammaSwap(Z), Qg = GammaSwap(Q);
            //@formatter:off
            R4 <<
               cdual(X.cwiseProduct(Xg).sum()), cdual(Y.cwiseProduct(Xg).sum()), cdual(Z.cwiseProduct(Xg).sum()), cdual(Q.cwiseProduct(Xg).sum()),
               cdual(X.cwiseProduct(Yg).sum()), cdual(Y.cwiseProduct(Yg).sum()), cdual(Z.cwiseProduct(Yg).sum()), cdual(Q.cwiseProduct(Yg).sum()),
               cdual(X.cwiseProduct(Zg).sum()), cdual(Y.cwiseProduct(Zg).sum()), cdual(Z.cwiseProduct(Zg).sum()), cdual(Q.cwiseProduct(Zg).sum()),
               cdual(X.cwiseProduct(Qg).sum()), cdual(Y.cwiseProduct(Qg).sum()), cdual(Z.cwiseProduct(Qg).sum()), cdual(Q.cwiseProduct(Qg).sum());
            //@formatter:on
        }
        Eigen::Matrix4cd r = de_dual(R4);

        Eigen::VectorXd params = Eigen::VectorXd::Ones(17);


        optim::algo_settings_t opt_set;
        opt_set.print_level = 1;
        bool success = optim::bfgs(params, opt_fn4, nullptr, opt_set);

        if (success) {
            std::cout << "bfgs: autodiff test completed successfully." << std::endl;
        } else {
            std::cout << "bfgs: autodiff test completed unsuccessfully." << std::endl;
        }

        Eigen::Matrix4cd A = Eigen::Matrix4cd::Zero();
        //@formatter:off
        A <<
          (i * params[0]),   (params[1] + i * params[2]),    (params[3] + i * params[4]), (params[5] + i * params[6]),
          (-params[1] + i * params[2]), (i * params[7]),     (params[8] + i * params[9]), (params[10] + i * params[11]),
          (-params[3] + i * params[4]), (-params[8] + i * params[9]),   (i * params[12]), (params[13] + i * params[14]),
          (-params[5] + i * params[6]), (-params[10] + i * params[11]), (-params[13] + i * params[14]), (i *params[15]);
        //@formatter:on
        A = std::exp(cd(0, params[16])) *
                ((A - Eigen::Matrix4cd::Identity()) * (Eigen::Matrix4cd::Identity() + A).inverse()).eval();
        std::cout << "Err_in = " << -(r).real().trace() << std::endl;
        std::cout << "Err_out = " << -(A.transpose() * r * A).real().trace() << std::endl;
        std::cout << "R  = " << r << std::endl;
        std::cout << "A^T * R * A = " << -A.transpose() * r * A << std::endl;
        std::cout << "A  = " << A << std::endl;
        //@formatter:off
        return std::make_tuple((A(0, 0) * X + A(0, 1) * Y + A(0, 2) * Z, A(0, 3) * Q).eval(),
                               (A(1, 0) * X + A(1, 1) * Y + A(1, 2) * Z, A(1, 3) * Q).eval(),
                               (A(2, 0) * X + A(2, 1) * Y + A(2, 2) * Z, A(2, 3) * Q).eval(),
                               (A(3, 0) * X + A(3, 1) * Y + A(3, 2) * Z, A(3, 3) * Q).eval());
        //@formatter:on
    }

/// Uses the Majoranize function to optimize Majorana Operators and then Builds according Fermi Operators
/// \param X Input Vector
/// \param Y Input Vector
/// \return tupel(A,B) Fermi Vectors
    std::tuple<Vec_cd, Vec_cd> Fermiize(const Vec_cd &X, const Vec_cd &Y) {
        auto tpl = Majoranize(X, Y);
        return std::make_tuple(1.0 / 2.0 * (std::get<0>(tpl) + cd(0, 1) * std::get<1>(tpl)),
                               1.0 / 2.0 * (std::get<0>(tpl) - cd(0, 1) * std::get<1>(tpl)));
    }

///
/// \param X
/// \param Y
/// \return
    std::tuple<Vec_cd, Vec_cd, Vec_cd, Vec_cd> FermiiMajize(const Vec_cd &X, const Vec_cd &Y) {
        auto tpl = Majoranize(X, Y);
        return std::make_tuple(std::get<0>(tpl), std::get<1>(tpl),
                               1.0 / 2.0 * (std::get<0>(tpl) + cd(0, 1) * std::get<1>(tpl)),
                               1.0 / 2.0 * (std::get<0>(tpl) - cd(0, 1) * std::get<1>(tpl)));
    }

    /// Uses the Majoranize function to optimize Majorana Operators and then Builds according Fermi Operators
/// \param X Input Vector
/// \param Y Input Vector
/// \return tupel(A,B) Fermi Vectors
    std::tuple<Vec_cd, Vec_cd, Vec_cd, Vec_cd>
    Fermiize(const Vec_cd &X, const Vec_cd &Y, const Vec_cd &Z, const Vec_cd &Q) {
        auto tpl = Majoranize(X, Y, Z, Q);
        //@formatter:off
        return std::make_tuple(1.0 / 2.0 * (std::get<0>(tpl) + cd(0, 1) * std::get<1>(tpl)),
                               1.0 / 2.0 * (std::get<0>(tpl) - cd(0, 1) * std::get<1>(tpl)),
                               1.0 / 2.0 * (std::get<2>(tpl) + cd(0, 1) * std::get<3>(tpl)),
                               1.0 / 2.0 * (std::get<2>(tpl) - cd(0, 1) * std::get<3>(tpl)));
        //@formatter:on
    }

///
/// \param X
/// \param Y
/// \return
    std::tuple<Vec_cd, Vec_cd, Vec_cd, Vec_cd, Vec_cd, Vec_cd, Vec_cd, Vec_cd, Vec_cd, Vec_cd, Vec_cd, Vec_cd, Vec_cd, Vec_cd>
    FermiiMajize(const Vec_cd &X, const Vec_cd &Y, const Vec_cd &Z, const Vec_cd &Q) {
        auto tpl = Majoranize(X, Y, Z, Q);
        //@formatter:off
        return std::make_tuple(std::get<0>(tpl), std::get<1>(tpl),
                               1.0 / 2.0 * (std::get<0>(tpl) + cd(0, 1) * std::get<1>(tpl)),
                               1.0 / 2.0 * (std::get<0>(tpl) - cd(0, 1) * std::get<1>(tpl)),
                               1.0 / 2.0 * (std::get<2>(tpl) + cd(0, 1) * std::get<3>(tpl)),
                               1.0 / 2.0 * (std::get<2>(tpl) - cd(0, 1) * std::get<3>(tpl)),
                               1.0 / 2.0 * (std::get<0>(tpl) + cd(0, 1) * std::get<3>(tpl)),
                               1.0 / 2.0 * (std::get<0>(tpl) - cd(0, 1) * std::get<3>(tpl)),
                               1.0 / 2.0 * (std::get<2>(tpl) + cd(0, 1) * std::get<1>(tpl)),
                               1.0 / 2.0 * (std::get<2>(tpl) - cd(0, 1) * std::get<1>(tpl)),
                               1.0 / 2.0 * (std::get<0>(tpl) + cd(0, 1) * std::get<2>(tpl)),
                               1.0 / 2.0 * (std::get<0>(tpl) - cd(0, 1) * std::get<2>(tpl)),
                               1.0 / 2.0 * (std::get<3>(tpl) + cd(0, 1) * std::get<1>(tpl)),
                               1.0 / 2.0 * (std::get<3>(tpl) - cd(0, 1) * std::get<1>(tpl)));
        //@formatter:on
    }

}