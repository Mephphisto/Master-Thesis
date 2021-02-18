#include <Eigen/Dense>
#include "Typedefs.hpp"

/// Error function of |gA^* -A|²+|gB^* -B|²
/// \param gxx <gX|X>
/// \param gyy <gY|Y>
/// \param re_gxy <gX|Y>
/// \param Phi Phi Angle
/// \param Xi Xi Angle
/// \param Zeta Zeata Angel
/// \param Theta Theta Angel
/// \return  |gA^* -A|²+|gB^* -B|²
double Error(cd gxx, cd gyy, double re_gxy, double Phi, double Xi, double Zeta, double Theta) {
    double CT = std::cos(Theta), C2T = CT * CT, ST = std::sin(Theta), S2T = ST * ST;
    return -gyy.real() * (C2T * std::cos(2.0 * Xi - 2.0 * Phi) +
                          std::cos(2.0 * Zeta - 2.0 * Phi) * S2T) -
           2.0 * gxx.real() * (C2T * std::cos(2.0 * (Zeta + Phi)) +
                               std::cos(2.0 * (Xi + Phi)) * S2T) +
           2.0 * gyy.imag() * (S2T * std::sin(2.0 * Zeta - 2.0 * Phi) +
                               C2T * std::sin(2.0 * Xi - 2.0 * Phi)) +
           2.0 * gxx.imag() * (C2T * std::sin(2.0 * (Zeta + Phi)) +
                               S2T * std::sin(2.0 * (Xi + Phi))) +
           2.0 * re_gxy * (std::cos(2.0 * Zeta) - std::cos(2.0 * Xi)) * std::sin(2.0 * Theta);
}

/// Computes teh Derivateive dError/dTheta
/// \param gxx <gX|X>
/// \param gyy <gY|Y>
/// \param re_gxy <gX|Y>
/// \param Phi Phi Angle
/// \param Xi Xi Angle
/// \param Zeta Zeata Angel
/// \param Theta Theta Angel
/// \return dError/dTheta
inline double dEdTheta(cd gxx, cd gyy, double re_gxy, double Phi, double Xi, double Zeta, double Theta) {
    return 4.0 * (std::cos(2.0 * Theta) * (std::cos(2.0 * Zeta) - std::cos(2.0 * Xi)) * re_gxy
                  + std::sin(2.0 * Theta) * (std::sin(Zeta - Xi) *
                                             (-std::cos(Zeta + Xi + 2.0 * Phi) * gxx.imag())
                                             + std::cos(Zeta + Xi - 2.0 * Phi) * gyy.imag()
                                             + std::sin(Zeta + Xi - 2.0 * Phi) * gyy.real()
                                             - std::sin(Zeta + Xi + 2.0 * Phi) * gxx.real()));
}

/// Computes teh Derivateive dError/dPhi
/// \param gxx <gX|X>
/// \param gyy <gY|Y>
/// \param re_gxy <gX|Y>
/// \param Phi Phi Angle
/// \param Xi Xi Angle
/// \param Zeta Zeata Angel
/// \param Theta Theta Angel
/// \return dError/dPhi
inline double dEdPhi(cd gxx, cd gyy, double re_gxy, double Phi, double Xi, double Zeta, double Theta) {
    double CT = std::cos(Theta), C2T = CT * CT, ST = std::sin(Theta), S2T = ST * ST;
    return 4.0 * (-gyy.imag() * (C2T * std::cos(2.0 * Xi - 2.0 * Phi) +
                                 std::cos(2.0 * Zeta - 2.0 * Phi) * S2T) +
                  gxx.imag() * (C2T * std::cos(2.0 * (Zeta + Phi)) +
                                std::cos(2.0 * (Xi + Phi)) * S2T) -
                  gyy.real() * S2T * std::sin(2.0 * Zeta - 2.0 * Phi) -
                  C2T * gyy.real() * std::sin(2.0 * Xi - 2.0 * Phi) +
                  C2T * gxx.real() * std::sin(2.0 * (Zeta + Phi)) +
                  gxx.real() * S2T * std::sin(2.0 * (Xi + Phi)));
}

/// Computes teh Derivateive dError/dXi
/// \param gxx <gX|X>
/// \param gyy <gY|Y>
/// \param re_gxy <gX|Y>
/// \param Phi Phi Angle
/// \param Xi Xi Angle
/// \param Zeta Zeata Angel
/// \param Theta Theta Angel
/// \return dError/dXI
inline double dEdXi(cd gxx, cd gyy, double re_gxy, double Phi, double Xi, double Zeta, double Theta) {
    double CT = std::cos(Theta), C2T = CT * CT, ST = std::sin(Theta), S2T = ST * ST;
    return 4.0 * (C2T * std::cos(2.0 * Xi - 2.0 * Phi) * gyy.imag() +
                  std::cos(2.0 * (Xi + Phi)) * gxx.imag() * S2T +
                  re_gxy * std::sin(2.0 * Theta) * std::sin(2.0 * Xi) +
                  C2T * gyy.real() * std::sin(2.0 * Xi - 2.0 * Phi) +
                  gxx.real() * S2T * std::sin(2.0 * (Xi + Phi)));
}

/// Computes teh Derivateive dError/dZeta
/// \param gxx <gX|X>
/// \param gyy <gY|Y>
/// \param re_gxy <gX|Y>
/// \param Phi Phi Angle
/// \param Xi Xi Angle
/// \param Zeta Zeata Angel
/// \param Theta Theta Angel
/// \return dError/dZeta
inline double dEdZeta(cd gxx, cd gyy, double re_gxy, double Phi, double Xi, double Zeta, double Theta) {
    double CT = std::cos(Theta), C2T = CT * CT, ST = std::sin(Theta), S2T = ST * ST;
    return 4.0 * (C2T * std::cos(2.0 * (Zeta + Phi)) * gxx.imag() +
                  std::cos(2.0 * Zeta - 2.0 * Phi) * gyy.imag() * S2T -
                  re_gxy * std::sin(2 * Zeta) * std::sin(2 * Theta) +
                  gyy.real() * S2T * std::sin(2.0 * Zeta - 2.0 * Phi) +
                  C2T * gxx.real() * std::sin(2.0 * (Zeta + Phi)));
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

/// Does GD update and maps back to [0, 2Pi]
/// \param X Value
/// \param dX derivative
/// \param lambda rearinig rate
/// \return new value
inline double upd_2Pi(double X, double dX, double lambda) {
    return std::fmod(X - lambda * dX, 2.0 * M_PI); // X - lambda * dX ;
    //;
}

/// This function Performs an GD optimisation to find superpositions of X and Y such that they fullfill (gamma A^* == A)
/// \param X Input Vector
/// \param Y Input Vector
/// \return tupel (A,B)
std::tuple<Vec_cd, Vec_cd> Majoranize(const Vec_cd &X, const Vec_cd &Y) {
    double lambda = std::numeric_limits<double>::max();
    double Phi = 0.0, Xi = 0.0, Zeta = 0.0, Theta = 0.0;
    assert(X.size() == Y.size());
    Vec_cd gX(GammaSwap(X)), gY(GammaSwap(Y));
#ifdef DEBUG_ACTIVE
    std::cout << "gX*X = " << gX.dot(X.conjugate()) << " gY*Y = " << gY.dot(Y.conjugate()) << std::endl;
#endif
    cd gxx = gX.conjugate().dot(X), gyy = gY.conjugate().dot(Y);
    double re_gxy = gX.conjugate().dot(Y).real();
    double dEdP, dEdX, dEdZ, dEdT;
#ifdef DEBUG_ACTIVE
    std::cout << "Beginnig Majoranize, Error_0 = " << 4 - 0 - 2.0 * Error(gxx, gyy, re_gxy, Phi, Xi, Zeta, Theta)
              << std::endl
              << " Teta = " << Theta << " Phi = " << Phi << " Xi = " << Xi << " Zeta = " << Zeta << std::endl
              << " Check 4-2Re(X*gX)-2Re(Y*gY) = " << X.dot(X) + Y.dot(Y) + gX.dot(gX) + gY.dot(gY) - 2.0 * (std::real(
            X.adjoint().dot(gX)) + std::real(Y.adjoint().dot(gY))) << std::endl
              << " Check |gX-X^*|^2 + |gY-Y^*|^2 = "
              << (gX - X.conjugate()).squaredNorm() + (gY - Y.conjugate()).squaredNorm() << std::endl
              << " gxx = " << gxx << " gyy = " << gyy << " re_gxy = " << re_gxy
              << std::endl;
#endif
    for (size_t k = 0; k < 1000; k++) {
        double Error_Old = Error(gxx, gyy, re_gxy, Phi, Xi, Zeta, Theta);
        dEdP = dEdPhi(gxx, gyy, re_gxy, Phi, Xi, Zeta, Theta);
        dEdX = dEdXi(gxx, gyy, re_gxy, Phi, Xi, Zeta, Theta);
        dEdZ = dEdZeta(gxx, gyy, re_gxy, Phi, Xi, Zeta, Theta);
        dEdT = dEdTheta(gxx, gyy, re_gxy, Phi, Xi, Zeta, Theta);
        while (Error_Old <= Error(gxx, gyy, re_gxy,
                                  upd_2Pi(Phi, dEdP, lambda),
                                  upd_2Pi(Xi, dEdX, lambda),
                                  upd_2Pi(Zeta, dEdZ, lambda),
                                  upd_2Pi(Theta, dEdT, lambda))) {
            lambda /= 2;
            if (lambda < std::numeric_limits<double>::epsilon()) {
                k = 1000;
#ifdef DEBUG_ACTIVE
                std::cout << "Reached Convergence" << std::endl;
#endif
                break;
            }
        }
#ifdef DEBUG_ACTIVE
        std::cout << "Error Improvement " << Error_Old - Error(gxx, gyy, re_gxy,
                                                               upd_2Pi(Phi, dEdP, lambda),
                                                               upd_2Pi(Xi, dEdX, lambda),
                                                               upd_2Pi(Zeta, dEdZ, lambda),
                                                               upd_2Pi(Theta, dEdT, lambda)) << std::endl;
#endif

        Phi = upd_2Pi(Phi, dEdP, lambda);
        Xi = upd_2Pi(Xi, dEdX, lambda);
        Zeta = upd_2Pi(Zeta, dEdZ, lambda);
        Theta = upd_2Pi(Theta, dEdT, lambda);
#ifdef DEBUG_ACTIVE
        if (k >= 999) {
            std::cout << "Aborted with lambda = " << lambda << std::endl;
        }
#endif
    }
    Vec_cd A = (std::exp(cd(0, Xi)) * (std::exp(cd(0, Phi)) * std::sin(Theta) * X +
                                       std::exp(cd(0, -Phi)) * std::cos(Theta) * Y)).normalized(),
            B = std::exp(cd(0, Zeta)) * (std::exp(cd(0, Phi)) * std::cos(Theta) * X -
                                         std::exp(cd(0, -Phi)) * std::sin(Theta) * Y);
#ifdef DEBUG_ACTIVE
    std::cout << "Ending Majoranize, Error_Final = " << 4.0 - 2.0 * Error(gxx, gyy, re_gxy, Phi, Xi, Zeta, Theta)
              << std::endl
              << " Teta = " << Theta << " Phi = " << Phi << " Xi = " << Xi << " Zeta = " << Zeta << std::endl
              << " Check |gA-A^*|^2 + |gB-B^*|^2 = "
              << (GammaSwap(A) - A.conjugate()).squaredNorm() + (GammaSwap(B) - B.conjugate()).squaredNorm()
              << std::endl;
#endif
    return
            std::make_tuple(A, B);
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