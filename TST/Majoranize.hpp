#include <Eigen/Dense>
#include "Typedefs.hpp"

/// Error function
/// \param gxx_gyy (gamma.X).X+(gamma.Y).Y (scalar Product not inner product)
/// \param re_gxy (gamma.X).Y (scalar Product not inner product)
/// \param Phi_A Phi_A Angle
/// \param Phi_B Phi_B Angle
/// \param Theta Theta Angle
/// \return Error/2-4
double Error(cd gxx_gyy, double re_gxy, double Phi_A, double Phi_B, double Theta) {
    double ST = std::sin(Theta), CT = std::cos(Theta), S2_T = ST * ST, C2_T = CT * CT;
    return -((std::cos(2 * Phi_A) * S2_T + std::cos(2 * Phi_B) * C2_T) * gxx_gyy.real()
             - (std::sin(2 * Phi_A) * S2_T + std::sin(2 * Phi_B) * C2_T) * gxx_gyy.imag()
             + 4.0 * std::cos(Phi_A + Phi_B) * CT * ST * re_gxy);
}

/// Derivatieve of Error function in respect of Theta
/// \param gxx_gyy (gamma.X).X+(gamma.Y).Y (scalar Product not inner product)
/// \param re_gxy (gamma.X).Y (scalar Product not inner product)
/// \param Phi_A Phi_A Angle
/// \param Phi_B Phi_B Angle
/// \param Theta Theta Angle
/// \return dError/dTheta
inline double dEdTheta(cd gxx_gyy, double re_gxy, double Phi_A, double Phi_B, double Theta) {
    double ST = std::sin(Theta), CT = std::cos(Theta), S2_T = ST * ST, C2_T = CT * CT;
    return 4.0 * (sin(2 * Theta) * std::sin(Phi_A - Phi_B) *
                  (
                          std::sin(Phi_A + Phi_B) * gxx_gyy.real()
                          + std::cos(Phi_A + Phi_B) * gxx_gyy.imag()
                  )
                  - 4.0 * std::cos(Phi_A + Phi_B) * std::cos(2 * Theta) * re_gxy);
}

/// Derivatieve of Error function in respect of Phi_A
/// \param gxx_gyy (gamma.X).X+(gamma.Y).Y (scalar Product not inner product)
/// \param re_gxy (gamma.X).Y (scalar Product not inner product)
/// \param Phi_A Phi_A Angle
/// \param Phi_B Phi_B Angle
/// \param Theta Theta Angle
/// \return dError/dPhi_A
inline double dEdPhi_A(cd gxx_gyy, double re_gxy, double Phi_A, double Phi_B, double Theta) {
    double ST = std::sin(Theta), CT = std::cos(Theta), S2_T = ST * ST;
    return 4.0 * ST * (std::sin(2 * Phi_A) * ST * gxx_gyy.real()
                       + std::cos(2 * Phi_A) * ST * gxx_gyy.imag()
                       + 2.0 * std::sin(Phi_A + Phi_B) * CT * re_gxy);
}

/// Derivatieve of Error function in respect of Phi_B
/// \param gxx_gyy (gamma.X).X+(gamma.Y).Y (scalar Product not inner product)
/// \param re_gxy (gamma.X).Y (scalar Product not inner product)
/// \param Phi_A Phi_A Angle
/// \param Phi_B Phi_B Angle
/// \param Theta Theta Angle
/// \return dError/dPhi_B
inline double dEdPhi_B(cd gxx_gyy, double re_gxy, double Phi_A, double Phi_B, double Theta) {
    double ST = std::sin(Theta), CT = std::cos(Theta), C2_T = CT * CT;
    return 4.0 * CT * (std::sin(2 * Phi_B) * CT * gxx_gyy.real()
                       + std::cos(2 * Phi_B) * CT * gxx_gyy.imag()
                       + 2.0 * std::sin(Phi_A + Phi_B) * ST * re_gxy);
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
    return std::fmod(X - lambda * dX, 2 * M_PI); // X - lambda * dX ;
    //;
}

/// This function Performs an GD optimisation to find superpositions of X and Y such that they fullfill (gamma A^* == A)
/// \param X Input Vector
/// \param Y Input Vector
/// \return tupel (A,B)
std::tuple<Vec_cd, Vec_cd> Majoranize(const Vec_cd &X, const Vec_cd &Y) {
    double lambda = 10000.0;
    double Phi_A = 0.0, Phi_B = 0.0, Theta = 0.0;
    assert(X.size() == Y.size());
    Vec_cd gX(GammaSwap(X)), gY(GammaSwap(Y));
#ifdef DEBUG_ACTIVE
    std::cout << "gX*X = " << gX.dot(X.conjugate()) << " gY*Y = " << gY.dot(Y.conjugate()) << std::endl;
#endif
    cd gxx_gyy = gX.conjugate().dot(X) + gY.conjugate().dot(Y);
    const size_t N = X.size();
    double re_gxy = gX.conjugate().dot(Y).real();
    /*for (size_t k = 0; k < X.size(); k++) {
        gxx_gyy += X[k] * X[(k + N) % (2 * N)] + Y[k] * Y[(k + N) % (2 * N)];
        re_gxy += std::real(X[(k + N) % (2 * N)] * Y[k]);
    }*/

    double dEdPA, dEdPB, dEdPT;
#ifdef DEBUG_ACTIVE
    std::cout << "Beginnig Majoranize, Error_0 = " << 4.0 + 2.0 * Error(gxx_gyy, re_gxy, Phi_A, Phi_B, Theta)
              << std::endl
              << " Teta = " << Theta << " Phi_A = " << Phi_A << " Phi_B = " << Phi_B << std::endl
              << " Check 4-2Re(X*gX)-2Re(Y*gY) = " << X.dot(X) + Y.dot(Y) + gX.dot(gX) + gY.dot(gY) - 2.0 * (std::real(
            X.adjoint().dot(gX)) + std::real(Y.adjoint().dot(gY))) << std::endl
              << " Check |gX-X^*|^2 + |gY-Y^*|^2 = "
              << (gX - X.conjugate()).squaredNorm() + (gY - Y.conjugate()).squaredNorm() << std::endl
              << " gxx_gyy = " << gxx_gyy << " re_gxy = " << re_gxy
              << std::endl;
#endif
    for (size_t k = 0; k < 100; k++) {
        double Error_Old = Error(gxx_gyy, re_gxy, Phi_A, Phi_B, Theta);
        dEdPA = dEdPhi_A(gxx_gyy, re_gxy, Phi_A, Phi_B, Theta);
        dEdPB = dEdPhi_B(gxx_gyy, re_gxy, Phi_A, Phi_B, Theta);
        dEdPT = dEdTheta(gxx_gyy, re_gxy, Phi_A, Phi_B, Theta);
        while (Error_Old < Error(gxx_gyy, re_gxy,
                                 upd_2Pi(Phi_A, dEdPA, lambda),
                                 upd_2Pi(Phi_B, dEdPB, lambda),
                                 upd_2Pi(Theta, dEdPT, lambda))) {
            lambda /= 2;
            if (lambda < std::numeric_limits<double>::epsilon()) {
                k = 200;
#ifdef DEBUG_ACTIVE
                std::cout << "Reached Convergence" << std::endl;
#endif
                break;
            }
        }
#ifdef DEBUG_ACTIVE
        std::cout << "Error Improvement " << Error_Old - Error(gxx_gyy, re_gxy,
                                                               upd_2Pi(Phi_A, dEdPA, lambda),
                                                               upd_2Pi(Phi_B, dEdPB, lambda),
                                                               upd_2Pi(Theta, dEdPT, lambda)) << std::endl;
#endif

        Phi_A = upd_2Pi(Phi_A, dEdPA, lambda);
        Phi_B = upd_2Pi(Phi_B, dEdPB, lambda);
        Theta = upd_2Pi(Theta, dEdPT, lambda);
#ifdef DEBUG_ACTIVE
        if (k == 99) {
            std::cout << "Aborted with lambda = " << lambda << std::endl;
        }
#endif
    }

    cd a = std::exp(cd(0, Phi_A)) * std::sin(Theta), b = std::exp(cd(0, Phi_B)) * std::cos(Theta);
    Vec_cd A = (a * X + b * Y).normalized(), B = (b * X + a * Y).normalized();
#ifdef DEBUG_ACTIVE
    std::cout << "Ending Majoranize, Error_Final = " << 4.0 + 2.0 * Error(gxx_gyy, re_gxy, Phi_A, Phi_B, Theta)
              << std::endl
              << " Teta = " << Theta << " Phi_A = " << Phi_A << " Phi_B = " << Phi_B << std::endl
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