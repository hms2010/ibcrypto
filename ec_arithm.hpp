#pragma once

#include "fp2_arithm.hpp"

template <const big_int_t& p> 
struct MontgomeryPoint
{
        // Montgomery curve point on projective space
        // Instead of (X:Y:Z) use (X:Z)
        Fp2Elem<p> X;
        Fp2Elem<p> Z;

        MontgomeryPoint(void); // Build T=(0:0:1)=(0:0) here
        ~MontgomeryPoint(void) = default;
        MontgomeryPoint(const Fp2Elem<p>& X);
        MontgomeryPoint(const Fp2Elem<p>& X, const Fp2Elem<p>& Z);
        MontgomeryPoint(const MontgomeryPoint& src);
        MontgomeryPoint(MontgomeryPoint&& src);

        MontgomeryPoint<p>& operator = (const MontgomeryPoint& src);
        MontgomeryPoint<p>& operator = (MontgomeryPoint&& src);

        MontgomeryPoint<p> get_normalized(void) const;
        Fp2Elem<p> get_x(void) const; // x = X / Z

        void print(void) const;
};

template <const big_int_t& p> 
class MontgomeryCurve
{
        // b*y^2 = x^3 + a*x^2 + x
        // a = A / C, b = B / C
        // B*y^2 = C * x^3 + A * x^2 +  C * x
        // x = X / Z, y = Y / Z
        // B*(Y^2)*Z = C*X^3 + A*(X^2)* Z + C*X*Z^2
        Fp2Elem<p> A;
        Fp2Elem<p> B;
        Fp2Elem<p> C;

        Fp2Elem<p> Ap24;
        Fp2Elem<p> Am24;
        Fp2Elem<p> a24;
        Fp2Elem<p> C24;
    public:
        MontgomeryCurve(void)=default;
        ~MontgomeryCurve(void)=default;

        MontgomeryCurve(const Fp2Elem<p>& A, const Fp2Elem<p>& B, const Fp2Elem<p>& C);
        MontgomeryCurve(const MontgomeryPoint<p>& P, const MontgomeryPoint<p>& Q, const MontgomeryPoint<p>& R);
        MontgomeryCurve(const MontgomeryCurve& src);
        MontgomeryCurve(MontgomeryCurve&& src);

        MontgomeryCurve<p>& operator = (const MontgomeryCurve<p>& src);
        MontgomeryCurve<p>& operator = (MontgomeryCurve<p>&& src);

        void xDbl(const MontgomeryPoint<p>& P, MontgomeryPoint<p>& R) const;
        void xTrpl(const MontgomeryPoint<p>& P, MontgomeryPoint<p>& R) const;
        void xDbl_e(const FpElem<p>& e, const MontgomeryPoint<p>& P, MontgomeryPoint<p>& R) const;
        void xTrpl_e(const FpElem<p>& e, const MontgomeryPoint<p>& P, MontgomeryPoint<p>& R) const;

        void xDblAdd(const MontgomeryPoint<p>& P, const MontgomeryPoint<p>& Q, const MontgomeryPoint<p>& R, 
                        MontgomeryPoint<p>& P2, MontgomeryPoint<p>& PQ) const;

        void ladder3pt(const FpElem<p> k, const MontgomeryPoint<p>& P, const MontgomeryPoint<p>& Q, const MontgomeryPoint<p>& R,
                        MontgomeryPoint<p>& S) const;

        void iso2_curve(const MontgomeryPoint<p>& G, MontgomeryCurve<p>& curve) const;
        void iso2_eval(const MontgomeryPoint<p>& G, const MontgomeryPoint<p>& P, MontgomeryPoint<p>& R) const;
        void iso2e(const FpElem<p>& e, const MontgomeryPoint<p>& G0, MontgomeryCurve<p>& curve, 
                        const MontgomeryPoint<p>& P, const MontgomeryPoint<p>& Q, const MontgomeryPoint<p>& R,
                              MontgomeryPoint<p>& P1,      MontgomeryPoint<p>& Q1,      MontgomeryPoint<p>& R1) const;

        void iso3_curve(const MontgomeryPoint<p>& G, MontgomeryCurve<p>& curve, Fp2Elem<p>& K1, Fp2Elem<p>& K2) const;
        void iso3_eval(const MontgomeryPoint<p>& G, const MontgomeryPoint<p>& P, const Fp2Elem<p>& K1, const Fp2Elem<p>& K2,
                        MontgomeryPoint<p>& R) const;
        void iso3e(const FpElem<p>& e, const MontgomeryPoint<p>& G0, MontgomeryCurve<p>& curve, 
                        const MontgomeryPoint<p>& P, const MontgomeryPoint<p>& Q, const MontgomeryPoint<p>& R,
                              MontgomeryPoint<p>& P1,      MontgomeryPoint<p>& Q1,      MontgomeryPoint<p>& R1) const;
        Fp2Elem<p> j_inv(void) const;
        Fp2Elem<p> get_A(void) const;
        Fp2Elem<p> get_C(void) const;
        Fp2Elem<p> get_a(void) const;

        void print(void) const;

};

#include "ec_arithm_impl.cpp"
