#include <Eigen/Dense>
#include "Typedefs.hpp"

double Error(cd gxx_gyy, double re_gxy, double Phi_A, double Phi_B, double Theta) {
    double ST = std::sin(Theta), CT = std::cos(Theta), S2_T = ST * ST, C2_T = CT * CT;
    return 4.0 - 2.0 * (std::real((std::exp(cd(0, 2 * Phi_A)) * S2_T + std::exp(cd(0, 2 * Phi_B)) * C2_T) * gxx_gyy)
                        + 2.0 * std::real(std::exp(cd(0, Phi_A + Phi_B)) * ST * CT * re_gxy));
}

inline double dEdTheta(cd gxx_gyy, double re_gxy, double Phi_A, double Phi_B, double Theta) {
    double ST = std::sin(Theta), CT = std::cos(Theta), S2_T = ST * ST, C2_T = CT * CT;
    return 2 * (ST * CT * (std::cos(2 * Phi_A) - std::cos(2 * Phi_B) * gxx_gyy.real() -
                           std::sin(2 * Phi_A) - std::sin(2 * Phi_B) * gxx_gyy.imag())
                + (C2_T - S2_T) * std::cos(Phi_A + Phi_B) * re_gxy);
}

inline double dEdPhi_A(cd gxx_gyy, double re_gxy, double Phi_A, double Phi_B, double Theta) {
    double ST = std::sin(Theta), CT = std::cos(Theta), S2_T = ST * ST;
    return -2.0 * (std::sin(2 * Phi_A) * S2_T * gxx_gyy.real()
                   + std::cos(2 * Phi_A) * S2_T * gxx_gyy.imag()
                   + std::sin(Phi_A + Phi_B) * ST * CT * re_gxy);
}

inline double dEdPhi_B(cd gxx_gyy, double re_gxy, double Phi_A, double Phi_B, double Theta) {
    double ST = std::sin(Theta), CT = std::cos(Theta), C2_T = CT * CT;
    return -2.0 * (std::sin(2 * Phi_B) * C2_T * gxx_gyy.real()
                   + std::cos(2 * Phi_B) * C2_T * gxx_gyy.imag()
                   + std::sin(Phi_A + Phi_B) * ST * CT * re_gxy);
}

Vec_cd GammaSwap(const Vec_cd &X) {
    const size_t N = X.size() / 2;
    Vec_cd gX(2 * N);
    gX.head(N) = X.tail(N);
    gX.tail(N) = X.head(N);
    return gX;
}

std::tuple<Vec_cd, Vec_cd> Majoranize(const Vec_cd &X, const Vec_cd &Y) {
    double lambda = 1.0;
    double Phi_A = 0.0, Phi_B = 0.0, Theta = M_PI / 2.0;
    assert(X.size() == Y.size());
    Vec_cd gX(GammaSwap(X)), gY(GammaSwap(Y));
#ifdef DEBUG_ACTIVE
    std::cout << "gX*X = " << gX.dot(X.conjugate()) << " gY*Y = " << gY.dot(Y.conjugate()) << std::endl;
#endif
    cd gxx_gyy = 0.0;
    const size_t N = X.size();
    double re_gxy = 0.0;
    for (size_t k = 0; k < X.size(); k++) {
        gxx_gyy += X[k] * X[(k + N) % (2 * N)] + Y[k] * Y[(k + N) % (2 * N)];
        re_gxy += std::real(X[(k + N) % (2 * N)] * Y[k]);
    }

    double dEdPA, dEdPB, dEdPT;
#ifdef DEBUG_ACTIVE
    std::cout << "Beginnig Majoranize, Error_0 = " << Error(gxx_gyy, re_gxy, 0.0, 0.0, 0.0) << std::endl
              << " Check 4-2Re(X*gX)-2Re(Y*gY) = " << X.dot(X) + Y.dot(Y) + gX.dot(gX) + gY.dot(gY) - 2.0 * (std::real(
            X.adjoint().dot(gX)) + std::real(Y.adjoint().dot(gY))) << std::endl
              << " Check |gX-X^*|^2 + |gY-Y^*|^2 "
              << (gX.conjugate() - X).squaredNorm() + (gY.conjugate() - Y).squaredNorm() << std::endl
              << " gxx_gyy = " << gxx_gyy << " re_gxy = " << re_gxy
              << std::endl;
#endif
    for (size_t k = 0; k < 100; k++) {
        double Error_Old = Error(gxx_gyy, re_gxy, Phi_A, Phi_B, Theta);
        dEdPA = dEdPhi_A(gxx_gyy, re_gxy, Phi_A, Phi_B, Theta);
        dEdPB = dEdPhi_B(gxx_gyy, re_gxy, Phi_A, Phi_B, Theta);
        dEdPT = dEdTheta(gxx_gyy, re_gxy, Phi_A, Phi_B, Theta);
        while (Error(gxx_gyy, re_gxy, Phi_A - lambda * dEdPA, Phi_B - lambda * dEdPB, Theta - lambda * dEdPT) >
               Error_Old) {
            lambda /= 2;
            if (lambda < std::numeric_limits<double>::epsilon()) {
                k = 200;
                break;
            }
        }
        Phi_A -= lambda * dEdPA;
        Phi_B -= lambda * dEdPB;
        Theta -= lambda * dEdPT;

    }

    cd a = std::exp(cd(0, Phi_A)) * std::sin(Theta), b = std::exp(cd(0, Phi_B)) * std::cos(Theta);
    Vec_cd A = (a * X + b * Y), B = (b * X + a * Y);
#ifdef DEBUG_ACTIVE
    std::cout << "Ending Majoranize, Error_Final = " << Error(gxx_gyy, re_gxy, Phi_A, Phi_B, Theta)
              << " Check |gX-X^*|^2 + |gY-Y^*|^2 "
              << (GammaSwap(A) - A.conjugate()).squaredNorm() + (GammaSwap(B) - B.conjugate()).squaredNorm()
              << std::endl;
#endif
    return
            std::make_tuple(A, B);
}

std::tuple<Vec_cd, Vec_cd> Fermiize(const Vec_cd &X, const Vec_cd &Y) {
    auto tpl = Majoranize(X, Y);
    return std::make_tuple(1.0 / 2.0 * (std::get<0>(tpl) + cd(0, 1) * std::get<1>(tpl)),
                           1.0 / 2.0 * (std::get<0>(tpl) - cd(0, 1) * std::get<1>(tpl)));
}