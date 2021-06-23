#pragma once

#include "fp2_arithm.hpp"
#include "ec_arithm.hpp"


struct TorsionGroup
{
    std::string e;
    std::string xP[2];
    std::string xQ[2];
    std::string xR[2];
};

struct CurveCoeffs
{
    std::string A[2];
    std::string B[2];
    std::string C[2];
};

enum ForsythiaSideParams
{
    for_side_2 = 0,
    for_side_3
};

enum _RealImage_
{
    real = 0,
    imag
};

enum ForsythiaSet
{
    forsythia80 = 0,
    forsythia128
};

struct ForsythiaParamSet
{
    std::string p;
    CurveCoeffs curve_coeffs;
    TorsionGroup torsion_params[2];
};

template <const big_int_t& p>
class Forsythia
{
        MontgomeryCurve<p> start_curve;

        ForsythiaSideParams side;
        FpElem<p> e;
        FpElem<p> sk;
        MontgomeryPoint<p> P;
        MontgomeryPoint<p> Q;
        MontgomeryPoint<p> R;

        MontgomeryPoint<p> nP;
        MontgomeryPoint<p> nQ;
        MontgomeryPoint<p> nR;

    public:
        Forsythia(void) = delete;
        ~Forsythia(void) = default;
        Forsythia(ForsythiaSet param_set, ForsythiaSideParams side);
        void isogen(const FpElem<p>& sk, MontgomeryPoint<p>& pkP, MontgomeryPoint<p>& pkQ, MontgomeryPoint<p>& pkR);
        void isoex(const FpElem<p>& sk, const MontgomeryPoint<p>& pkP, const MontgomeryPoint<p>& pkQ, const MontgomeryPoint<p>& pkR, Fp2Elem<p>& j_inv);
};






extern ForsythiaParamSet g_param_set[2];

template <const mpz_class& p>
Forsythia<p>::Forsythia(ForsythiaSet param_set, ForsythiaSideParams _side_)
{
    side = _side_;
    start_curve = MontgomeryCurve<p>(Fp2Elem<p>(g_param_set[param_set].curve_coeffs.A[real],
                                                g_param_set[param_set].curve_coeffs.A[imag]),
                                    Fp2Elem<p>(g_param_set[param_set].curve_coeffs.B[real],
                                               g_param_set[param_set].curve_coeffs.B[imag]),
                                    Fp2Elem<p>(g_param_set[param_set].curve_coeffs.C[real],
                                               g_param_set[param_set].curve_coeffs.C[imag]));

    e = FpElem<p>(g_param_set[param_set].torsion_params[side].e);
    P = MontgomeryPoint<p>(Fp2Elem<p>(g_param_set[param_set].torsion_params[side].xP[real],
                                      g_param_set[param_set].torsion_params[side].xP[imag]));
    Q = MontgomeryPoint<p>(Fp2Elem<p>(g_param_set[param_set].torsion_params[side].xQ[real],
                                      g_param_set[param_set].torsion_params[side].xQ[imag]));
    R = MontgomeryPoint<p>(Fp2Elem<p>(g_param_set[param_set].torsion_params[side].xR[real],
                                      g_param_set[param_set].torsion_params[side].xR[imag]));

    nP = MontgomeryPoint<p>(Fp2Elem<p>(g_param_set[param_set].torsion_params[1 - static_cast<int>(side)].xP[real],
                                       g_param_set[param_set].torsion_params[1 - static_cast<int>(side)].xP[imag]));
    nQ = MontgomeryPoint<p>(Fp2Elem<p>(g_param_set[param_set].torsion_params[1 - static_cast<int>(side)].xQ[real],
                                       g_param_set[param_set].torsion_params[1 - static_cast<int>(side)].xQ[imag]));
    nR = MontgomeryPoint<p>(Fp2Elem<p>(g_param_set[param_set].torsion_params[1 - static_cast<int>(side)].xR[real],
                                       g_param_set[param_set].torsion_params[1 - static_cast<int>(side)].xR[imag]));
    #ifdef DEBUG
    std::cout << "P ="; P.get_x().print();
    std::cout << "Q ="; Q.get_x().print();
    std::cout << "R ="; R.get_x().print();

    std::cout << "nP ="; nP.get_x().print();
    std::cout << "nQ ="; nQ.get_x().print();
    std::cout << "nR ="; nR.get_x().print();
    #endif // DEBUG
}

template <const mpz_class& p>
void Forsythia<p>::isogen(const FpElem<p>& sk, MontgomeryPoint<p>& pkP, MontgomeryPoint<p>& pkQ, MontgomeryPoint<p>& pkR)
{
    MontgomeryPoint<p> G;
    MontgomeryCurve<p> curve;

    start_curve.ladder3pt(sk, P, Q, R, G);

    #ifdef DEBUG
    std::cout << "G = "; G.get_x().print();
    #endif

    if(side == for_side_2)
        start_curve.iso2e(e, G, curve, nP, nQ, nR, pkP, pkQ, pkR);
    else // for_side_3
        start_curve.iso3e(e, G, curve, nP, nQ, nR, pkP, pkQ, pkR);

    #ifdef DEBUG
    curve.print();
    std::cout << "pkP ="; pkP.get_x().print();
    std::cout << "pkQ ="; pkQ.get_x().print();
    std::cout << "pkR ="; pkR.get_x().print();
    #endif // DEBUG
}

template <const mpz_class& p>
void Forsythia<p>::isoex(const FpElem<p>& sk, const MontgomeryPoint<p>& pkP, const MontgomeryPoint<p>& pkQ, const MontgomeryPoint<p>& pkR, Fp2Elem<p>& j_inv)
{
    MontgomeryPoint<p> G;
    MontgomeryCurve<p> curve(pkP, pkQ, pkR);
    MontgomeryCurve<p> common_curve;
    MontgomeryPoint<p> plug;

    #ifdef DEBUG
    std::cout << "side" << 1 - static_cast<int>(side) << ": ";
    curve.print();
    #endif // DEBUG

    curve.ladder3pt(sk, pkP, pkQ, pkR, G);

    #ifdef DEBUG
    std::cout << "G = "; G.get_x().print();
    #endif // DEBUG

    if(side == for_side_2)
        curve.iso2e(e, G, common_curve, pkP, pkQ, pkR, plug, plug, plug);
    else // for_side_3
        curve.iso3e(e, G, common_curve, pkP, pkQ, pkR, plug, plug, plug);

    j_inv = common_curve.j_inv();

    #ifdef DEBUG
    common_curve.print();
    j_inv.print();
    #endif // DEBUG

}
