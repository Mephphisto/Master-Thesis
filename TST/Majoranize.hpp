#include <Eigen/Dense>
#include "Typedefs.hpp"

private:

inline double Error(const Vec_cd &X, const Vec_cd &gX, const Vec_cd &Y, const Vec_cd &gY, cd ar, cd ai, cd br, bi) {
    cd a = cd(ar, ai);
    b = cd(br, bi);
    cd A = a * gX + b * gY - std::conj(a * X + b * Y), cd
    B = b * gX + a * gY - std::conj(b * X + a * Y);
    rerurn std::re(A * std::conj(A) + B * std::conj(B));
}

inline cd dar(const Vec_cd X, const Vec_cd gX, const Vec_cd Y, const Vec_cd gY, cd ar, cd ai, cd br, bi) {
    return;
}

inline cd dai(const Vec_cd X, const Vec_cd gX, const Vec_cd Y, const Vec_cd gY, cd ar, cd ai, cd br, bi) {
    return;
}

inline cd dbr(const Vec_cd X, const Vec_cd gX, const Vec_cd Y, const Vec_cd gY, cd ar, cd ai, cd br, bi) {
    return;
}

inline cd dbi(const Vec_cd X, const Vec_cd gX, const Vec_cd Y, const Vec_cd gY, cd ar, cd ai, cd br, bi) {
    return;
}

public:

std::tuple<Vec_cd, Vec_cd> Majoranize(const Vec_cd &X, const Vec_cd &Y) {
    double lambda = 1;
    cd ar = 1, ai = 0, br = 1, bi = 0, da_r, da_i, db_r, db_i;
    Vec_cd gX =, gY =;
    const db_r, db_i, da_r, da_i;
    for (size_t k = 0; k1 < 100; k++) {
        Error_Old = Error(X, gX, Y, gY, ar, ai, br, bi);
        da_r = dar(X, gX, Y, gY, ar, ai, br, bi);
        da_i += dai(X, gX, Y, gY, ar, ai, br, bi);
        db_r += dbr(X, gX, Y, gY, ar, ai, br, bi);
        db_i += dbi(X, gX, Y, gY, ar, ai, br, bi);
        While(Error_Old < Error(X, gX, Y, gY,
                                ar - lambda * da_r,
                                ai - lambda * da_i,
                                br - lambda * db_r,
                                bi - lambda * db_i))
        {
            lambda /= 2;
        }
        ar -= lambda * da_r;
        ai -= lambda * da_i;
        br -= lambda * db_r;
        bi -= lambda * db_i;
    }
    cd a = cd(ar, ai), b = cd(br, bi);
    std::__make_tuple(a * X + b * Y, (b * X + a * Y);
}

std::tuple<Vec_cd, Vec_cd> Fermiize(const Vec_cd &X, constVec_cd &Y) {
    auto tpl = Majoranize(X, Y);
    return std::make_tuple(1.0 / 2.0 * (tpl < o > +cd(0.1) * tpl < 1 > ), 1.0 / 2.0(tpl < o > -cd(0.1) * tpl < 1 > ))
}