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

    Eigen::Matrix2cd dedual(const Eigen::Matrix2<cdual> &x) {
        Eigen::Matrix2cd y;
        for (int k = 0; k < 2; ++k) {
            for (int j = 0; j < 2; ++j) {
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
        R.normalize();

        Eigen::Matrix2cd r = dedual(R);
        Eigen::VectorXd params = Eigen::VectorXd::Ones(5);


        optim::algo_settings_t opt_set;
        opt_set.print_level = 1;
        bool success = optim::bfgs(params, opt_fn, nullptr, opt_set);

        if (success) {
            std::cout << "bfgs: autodiff test completed successfully." << std::endl;
        } else {
            std::cout << "bfgs: autodiff test completed unsuccessfully." << std::endl;
        }

        Eigen::Matrix2cd A = Eigen::Matrix2cd::Zero();
        A << i * params[0], params[1] + i * params[2], -params[1] + i * params[2], i * params[3];
        A = std::exp(cd(0, params[4])) *
            ((A - Eigen::Matrix2cd::Identity()) * (Eigen::Matrix2cd::Identity() + A).inverse()).eval();
        std::cout << "Err_in = " << -(r).real().trace() << std::endl;
        std::cout << "Err_out = " << -(A.transpose() * r * A).real().trace() << std::endl;
        std::cout << "R  = " << r << std::endl;
        std::cout << "A^T * R * A = " << -A.transpose() * r * A << std::endl;
        std::cout << "A  = " << A << std::endl;
        return std::make_tuple(A(0, 0) * X + A(0, 1) * Y, A(1, 0) * X + A(1, 1) * Y);
    }

    std::tuple<Vec_cd, Vec_cd, Vec_cd, Vec_cd>
    Majoranize(const Vec_cd &X, const Vec_cd &Y, const Vec_cd &Z, const Vec_cd &Q) {
        Id4 <<
            cdual(1), cdual(0), cdual(0), cdual(0)
        cdual(0), cdual(1), cdual(0), cdual(0)
        cdual(0), cdual(0), cdual(1), cdual(0)
        cdual(0), cdual(0), cdual(0), cdual(1);
        {
            Vec_cd Xg = GammaSwap(X), Yg = GammaSwap(Y), Zg = GammaSwap(Z), Qg = GammaSwap(Q);
            R44 <<
                cdual(X.cwiseProduct(Xg).sum()), cdual(Y.cwiseProduct(Xg).sum()), cdual(
                    Z.cwiseProduct(Xg).sum()), cdual(Q.cwiseProduct(Xg).sum()),
                    cdual(X.cwiseProduct(Yg).sum()),
                    cdual(Y.cwiseProduct(Yg).sum())
            cdual(Y.cwiseProduct(Qg).sum()), cdual(Y.cwiseProduct(Qg).sum()), cdual(Y.cwiseProduct(Qg).sum()), cdual(
                    Y.cwiseProduct(Qg).sum());
        }
        R.normalize();

        Eigen::Matrix2cd r = dedual(R);
        Eigen::VectorXd params = Eigen::VectorXd::Ones(5);


        optim::algo_settings_t opt_set;
        opt_set.print_level = 1;
        bool success = optim::bfgs(params, opt_fn, nullptr, opt_set);

        if (success) {
            std::cout << "bfgs: autodiff test completed successfully." << std::endl;
        } else {
            std::cout << "bfgs: autodiff test completed unsuccessfully." << std::endl;
        }

        Eigen::Matrix2cd A = Eigen::Matrix2cd::Zero();
        A << i * params[0], params[1] + i * params[2], -params[1] + i * params[2], i * params[3];
        A = std::exp(cd(0, params[4])) *
            ((A - Eigen::Matrix2cd::Identity()) * (Eigen::Matrix2cd::Identity() + A).inverse()).eval();
        std::cout << "Err_in = " << -(r).real().trace() << std::endl;
        std::cout << "Err_out = " << -(A.transpose() * r * A).real().trace() << std::endl;
        std::cout << "R  = " << r << std::endl;
        std::cout << "A^T * R * A = " << -A.transpose() * r * A << std::endl;
        std::cout << "A  = " << A << std::endl;
        return std::make_tuple(A(0, 0) * X + A(0, 1) * Y, A(1, 0) * X + A(1, 1) * Y);
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

    std::tuple<Vec_cd, Vec_cd, Vec_cd, Vec_cd> FermiiMajize(const Vec_cd &X, const Vec_cd &Y) {
        auto tpl = Majoranize(X, Y);
        return std::make_tuple(std::get<0>(tpl), std::get<1>(tpl),
                               1.0 / 2.0 * (std::get<0>(tpl) + cd(0, 1) * std::get<1>(tpl)),
                               1.0 / 2.0 * (std::get<0>(tpl) - cd(0, 1) * std::get<1>(tpl)));
    }
}