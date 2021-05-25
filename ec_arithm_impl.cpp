#include "ec_arithm.hpp"

#include <iostream>

#include <gmp.h>
#include <gmpxx.h>

#include "fp2_arithm.hpp"

template <const mpz_class& p>
MontgomeryPoint<p>::MontgomeryPoint(void) : X(Fp2Elem<p>("1", "0")), Z(Fp2Elem<p>("0", "0"))
{}

template <const mpz_class& p>
MontgomeryPoint<p>::MontgomeryPoint(const Fp2Elem<p>& _X_) : X(_X_), Z(Fp2Elem<p>("1", "0"))
{}

template <const mpz_class& p>
MontgomeryPoint<p>::MontgomeryPoint(const Fp2Elem<p>& _X_, const Fp2Elem<p>& _Z_) : X(_X_), Z(_Z_)
{}

template <const mpz_class& p>
MontgomeryPoint<p>::MontgomeryPoint(const MontgomeryPoint<p>& src) : X(src.X), Z(src.Z)
{}

template <const mpz_class& p>
MontgomeryPoint<p>::MontgomeryPoint(MontgomeryPoint<p>&& src) : X(std::move(src.X)), Z(std::move(src.Z))
{}

template <const mpz_class& p>
MontgomeryPoint<p>& MontgomeryPoint<p>::operator = (const MontgomeryPoint& src)
{
    X = src.X;
    Z = src.Z;
    return *this;
}

template <const mpz_class& p>
MontgomeryPoint<p>& MontgomeryPoint<p>::operator = (MontgomeryPoint&& src)
{
    X = std::move(src.X);
    Z = std::move(src.Z);
    return *this;
}

template <const mpz_class& p>
MontgomeryPoint<p> MontgomeryPoint<p>::get_normalized(void) const
{
    return MontgomeryPoint<p>(X * Z.invert(), Fp2Elem<p>("1", "0"));
}

template <const mpz_class& p>
Fp2Elem<p> MontgomeryPoint<p>::get_x(void) const // x = X / Z
{
    return X * Z.invert();
}

template <const mpz_class& p>
void MontgomeryPoint<p>::print(void) const
{
    std::cout << "(" << X.get_real_str() << " + i * " << X.get_imag_str() << ", " 
                << Z.get_real_str() << " + i * " << Z.get_imag_str() << ") mod "<< p << std::endl;
}

template <const mpz_class& p>
MontgomeryCurve<p>::MontgomeryCurve(const Fp2Elem<p>& _A_, const Fp2Elem<p>& _B_, const Fp2Elem<p>& _C_) : A(_A_), Ap24(_A_), Am24(_A_), B(_B_), C(_C_)
{
    C24 = C + C;               // 2C
    Ap24 += C24;               // A + 2C
    Am24 -= C24;               // A - 2C
    C24 += C24;                // 4C
    // a24 = Ap24 * C.invert();   // (A + 2C)/C = A/C + 2
    Fp2Elem<p> t("2", "0");
    a24 = A + t;
    t = (t + t).invert();
    a24 = a24 * t;
}

template <const mpz_class& p>
MontgomeryCurve<p>& MontgomeryCurve<p>::operator = (const MontgomeryCurve<p>& src)
{
    A = src.A;
    B = src.B;
    C = src.C;
    Ap24 = src.Ap24;
    Am24 = src.Am24;
    a24 = src.a24;
    C24 = src.C24;
    return *this;
}

template <const mpz_class& p>
MontgomeryCurve<p>& MontgomeryCurve<p>::operator = (MontgomeryCurve<p>&& src)
{
    A = std::move(src.A);
    B = std::move(src.B);
    C = std::move(src.C);
    a24 = std::move(src.a24);
    C24 = std::move(src.C24);
    Ap24 = std::move(src.Ap24);
    Am24 = std::move(src.Am24);
    return *this;
}

template <const mpz_class& p>
void MontgomeryCurve<p>::xDbl(const MontgomeryPoint<p>& P, MontgomeryPoint<p>& R) const
{
    Fp2Elem<p> t0 = P.X - P.Z; // XP - ZP
    Fp2Elem<p> t1 = P.X + P.Z; // XP + ZP
    t0 =t0 * t0;               // (XP - ZP)^2
    t1 =t1 * t1;               // (XP + ZP)^2
    R.Z = C24 * t0;            // 4C * (XP - ZP)^2
    R.X = R.Z * t1;            // 4C * (XP - ZP)^2 * (XP + ZP)^2
    t1 = t1 - t0;              // (XP + ZP)^2 - (XP - ZP)^2
    t0 = Ap24 * t1;            // (A + 2С) * ((XP + ZP)^2 - (XP - ZP)^2)
    R.Z = R.Z + t0;            // 4C * (XP - ZP)^2 + (A + 2С) * ((XP + ZP)^2 - (XP - ZP)^2)
    R.Z = R.Z * t1;            // (4C * (XP - ZP)^2 + (A + 2С) * ((XP + ZP)^2 - (XP - ZP)^2)) * ((XP + ZP)^2 - (XP - ZP)^2)
}

template <const mpz_class& p>
void MontgomeryCurve<p>::xTrpl(const MontgomeryPoint<p>& P, MontgomeryPoint<p>& R) const
{
    Fp2Elem<p> t0(P.X - P.Z); // XP - ZP
    Fp2Elem<p> t1(P.X + P.Z); // XP + ZP
    Fp2Elem<p> t2(t0 * t0);   // (XP - ZP)^2
    Fp2Elem<p> t3(t1 * t1);   // (XP + ZP)^2
    Fp2Elem<p> t4(t1 + t0);   // (XP + ZP) + (XP - ZP) = 2XP
    t0 = t1 - t0;             // 2ZP
    t1 = t4 * t4;             // 4(XP)^2
    t1 = t1 - t3;             // 4(XP)^2 - (XP + ZP)^2
    t1 = t1 - t2;             // 4(XP)^2 - 2 * (XP^2 + ZP^2) = 2 * (XP^2 - ZP^2)
    Fp2Elem<p> t5(t3 * Ap24); // (A + 2C) * (XP + ZP)^2
    t3 = t3 * t5;             // (A + 2C) * (XP + ZP)^4
    Fp2Elem<p> t6(t2 * Am24); // (A - 2C) * (XP - ZP)^2
    t2 = t2 * t6;             // (A - 2C) * (XP - ZP)^4
    t3 = t2 - t3;             // (A - 2C) * (XP - ZP)^4 - (A + 2C) * (XP + ZP)^4
    t2 = t5 - t6;             // (A + 2C) * (XP + ZP)^2 - (A - 2C) * (XP - ZP)^2
    t1 = t1 * t2;             // 
    t2 = t3 + t1;             //
    t2 = t2 * t2;             //
    R.X = t2 * t4;            //
    t1 = t3 - t1;             //
    t1 = t1 * t1;             //
    R.Z = t1 * t0;            //
}

template <const mpz_class& p>
void MontgomeryCurve<p>::xDbl_e(const FpElem<p>& e, const MontgomeryPoint<p>& P, MontgomeryPoint<p>& R) const
{
    MontgomeryPoint<p> rR;
    R = P;
    for (auto i = FpElem<p>("0"); i < e; ++i)
    {
        xDbl(R, rR);
        R = rR;
    }
}

template <const mpz_class& p>
void MontgomeryCurve<p>::xTrpl_e(const FpElem<p>& e, const MontgomeryPoint<p>& P, MontgomeryPoint<p>& R) const
{
    MontgomeryPoint<p> rR;
    R = P;
    for (auto i = FpElem<p>("0"); i < e; ++i)
    {
        xTrpl(R, rR);
        R = rR;
    }
}

template <const mpz_class& p>
void MontgomeryCurve<p>::xDblAdd(const MontgomeryPoint<p>& P, const MontgomeryPoint<p>& Q, const MontgomeryPoint<p>& R, MontgomeryPoint<p>& P2, MontgomeryPoint<p>& PQ) const
{
    Fp2Elem<p> t0(P.X + P.Z);  // XP + ZP
    Fp2Elem<p> t1(P.X - P.Z);  // XP - ZP
    P2.X = t0 * t0;            // (XP + ZP)^2
    P2.Z = t1 * t1;            // (XP - ZP)^2
    PQ.X  = Q.X + Q.Z;         // XQ + ZQ
    Fp2Elem<p> t2(Q.X - Q.Z);  // XQ - ZQ
    t0 = t0 * t2;              // (XP + ZP) * (XQ - ZQ)
    t1 = t1 * PQ.X;            // (XP - ZP) * (XQ + ZQ)
    t2 = P2.X - P2.Z;          // (XP + ZP)^2 - (XP - ZP)^2 = 4*XP*ZP
    P2.X = P2.X * P2.Z;        // (XP - ZP)^2 * (XP + ZP)^2
    PQ.X = a24 * t2;           // (A/C + 2) * 4*XP*ZP
    PQ.Z = t0 - t1;            // (XP + ZP) * (XQ - ZQ) - (XP - ZP) * (XQ + ZQ)
    P2.Z = P2.Z + PQ.X;        // (XP - ZP)^2 + (A/C + 2) * 4*XP*ZP
    PQ.X = t0 + t1;            // (XP + ZP) * (XQ - ZQ) + (XP - ZP) * (XQ + ZQ)
    P2.Z = P2.Z * t2;          // 4*XP*ZP * ((XP - ZP)^2 + (A/C + 2) * 4*XP*ZP)
    PQ.Z = PQ.Z * PQ.Z;        // ((XP + ZP) * (XQ - ZQ) - (XP - ZP) * (XQ + ZQ))^2
    PQ.X = PQ.X * PQ.X;        // ((XP + ZP) * (XQ - ZQ) + (XP - ZP) * (XQ + ZQ))^2
    PQ.Z = PQ.Z * R.X;         // XR * ((XP + ZP) * (XQ - ZQ) - (XP - ZP) * (XQ + ZQ))^2
    PQ.X = PQ.X * R.Z;         // ZR * ((XP + ZP) * (XQ - ZQ) + (XP - ZP) * (XQ + ZQ))^2
}

template <const mpz_class& p>
void MontgomeryCurve<p>::ladder3pt(const FpElem<p> k, const MontgomeryPoint<p>& P, const MontgomeryPoint<p>& Q, const MontgomeryPoint<p>& R, MontgomeryPoint<p>& S) const
{
    MontgomeryPoint<p> S0(Q.get_x());
    MontgomeryPoint<p> S1(P.get_x());
    MontgomeryPoint<p> S2(R.get_x());

    MontgomeryPoint<p> R0;
    MontgomeryPoint<p> R1;
    MontgomeryPoint<p> R2;

    FpElem<p> i = k;
    const FpElem<p> fp_null("0");
    const FpElem<p> fp_two("2");
    while (i > fp_null)
    {
        if(i % fp_two == fp_null)
        {
            xDblAdd(S0, S2, S1, R0, R2);
            S0 = R0;
            S2 = R2;
        }
        else
        {
            xDblAdd(S0, S1, S2, R0, R1);
            S0 = R0;
            S1 = R1;
        }
        i /= fp_two;
    }
    S = std::move(S1);
}

template <const mpz_class& p>
void MontgomeryCurve<p>::iso2_curve(const MontgomeryPoint<p>& G, MontgomeryCurve<p>& curve) const
{
    Fp2Elem<p> rA(G.X * G.X); // XG^2
    Fp2Elem<p> rC(G.Z * G.Z); // ZG^2
    rA = rA + rA;             // 2*XG^2
    rA = rC - rA;             // ZG^2 - 2*XG^2
    rA = rA + rA;             // 2(ZG^2 - 2*XG^2)
    curve = MontgomeryCurve<p>(rA, Fp2Elem<p>("1", "0"), rC);
}

template <const mpz_class& p>
void MontgomeryCurve<p>::iso2_eval(const MontgomeryPoint<p>& G, const MontgomeryPoint<p>& P, MontgomeryPoint<p>& R) const
{
    Fp2Elem<p> t0(G.X + G.Z); // XG + ZG
    Fp2Elem<p> t1(G.X - G.Z); // XG - ZG
    R.X = P.X + P.Z;          // XP + ZP
    R.Z = P.X - P.Z;          // XP - ZP
    t0 = t0 * R.Z;            // (XG + ZG) * (XP - ZP)
    t1 = t1 * R.X;            // (XG - ZG) * (XP + ZP)
    R.X = t0 + t1;            // (XG + ZG) * (XP - ZP) + (XG - ZG) * (XP + ZP)
    R.Z = t0 - t1;            // (XG + ZG) * (XP - ZP) - (XG - ZG) * (XP + ZP)
    R.X = R.X * P.X;          // XP * ((XG + ZG) * (XP - ZP) + (XG - ZG) * (XP + ZP))
    R.Z = R.Z * P.Z;          // ZP * ((XG + ZG) * (XP - ZP) - (XG - ZG) * (XP + ZP))
}

template <const mpz_class& p>
void MontgomeryCurve<p>::iso2e(const FpElem<p>& e, const MontgomeryPoint<p>& G0, MontgomeryCurve<p>& curve, 
                        const MontgomeryPoint<p>& P, const MontgomeryPoint<p>& Q, const MontgomeryPoint<p>& R,
                              MontgomeryPoint<p>& P1,      MontgomeryPoint<p>& Q1,      MontgomeryPoint<p>& R1) const
{
    MontgomeryPoint<p> G(G0.get_normalized());
    MontgomeryPoint<p> rG;
    MontgomeryPoint<p> rP1;
    MontgomeryPoint<p> rQ1;
    MontgomeryPoint<p> rR1;
    MontgomeryPoint<p> T;
    P1 = P.get_normalized();
    Q1 = Q.get_normalized();
    R1 = R.get_normalized();
    curve = *this;
    for (auto i = e - FpElem<p>("1"); i >= FpElem<p>("0"); --i)
    {
        curve.xDbl_e(i, G, T);       // Ti = [2^i]Gi
        curve.iso2_curve(T, curve);  // Ei+1 = Ei/<Ti>
        curve.iso2_eval(T, P1, rP1); // Pi+1= phi2_i(Pi)
        curve.iso2_eval(T, Q1, rQ1); // Qi+1= phi2_i(Qi)
        curve.iso2_eval(T, R1, rR1); // Ri+1= phi2_i(Ri)
        P1 = rP1;
        Q1 = rQ1;
        R1 = rR1;
        if(i == FpElem<p>("0"))
            return;
        curve.iso2_eval(T, G, rG);  // Gi+1 = phi2_i(Gi)
        G = rG;
    }
}

template <const mpz_class& p>
void MontgomeryCurve<p>::iso3_curve(const MontgomeryPoint<p>& G, MontgomeryCurve<p>& curve, Fp2Elem<p>& K1, Fp2Elem<p>& K2) const
{
    K1 = G.X - G.Z;         // XG - ZG
    K2 = G.X + G.Z;         // XG + ZG
    Fp2Elem<p> t0(K1 * K1); // K1^2
    Fp2Elem<p> t1(K2 * K2); // K2^2
    Fp2Elem<p> t2(t0 + t1); // K1^2 + K2^2
    Fp2Elem<p> t3(K1 + K2); // K1 + K2
    t3 = t3 * t3;           // (K1 + K2)^2
    t3 =t3 - t2;            // (K1 + K2)^2 - (K1^2 + K2^2)
    t2 = t1 + t3;           // K2^2 + (K1 + K2)^2 - (K1^2 + K2^2)
    t3 = t3 + t0;           // K1^2 + (K1 + K2)^2 - (K1^2 + K2^2)
    Fp2Elem<p> t4(t3 + t0); // 2*K1^2 + (K1 + K2)^2 - (K1^2 + K2^2)
    t4 = t4 + t4;           // 2 * (2*K1^2 + (K1 + K2)^2 - (K1^2 + K2^2))
    t4 = t4 + t1;           // 
    curve.Am24 = t2 * t4;   // 
    t4 = t1 + t2;
    t4 = t4 + t4;
    t4 = t4 + t0;
    curve.Ap24 = t3 * t4;
    curve.A = curve.Ap24 + curve.Ap24 + curve.Am24 + curve.Am24;
    curve.C = curve.Ap24 - curve.Am24;
    curve = MontgomeryCurve<p>(curve.A, Fp2Elem<p>("1", "0"), curve.C);
}

template <const mpz_class& p>
void MontgomeryCurve<p>::iso3_eval(const MontgomeryPoint<p>& G, const MontgomeryPoint<p>& P, const Fp2Elem<p>& K1, const Fp2Elem<p>& K2, MontgomeryPoint<p>& R) const
{
    Fp2Elem<p> t0(P.X + P.Z); // XP + ZP
    Fp2Elem<p> t1(P.X - P.Z); // XP - ZP
    t0 = t0 * K1;             // K1 * (XP + ZP)
    t1 = t1 * K2;             // K2 * (XP - ZP)
    Fp2Elem<p> t2(t0 + t1);   // K1 * (XP + ZP) + K2 * (XP - ZP)
    t0 = t1 - t0;             // K2 * (XP - ZP) - K1 * (XP + ZP)
    t2 = t2 * t2;             // (K1 * (XP + ZP) + K2 * (XP - ZP))^2
    t0 = t0 * t0;             // (K1 * (XP + ZP) - K2 * (XP - ZP))^2
    R.X = P.X * t2;           // XP * (K1 * (XP + ZP) + K2 * (XP - ZP))^2
    R.Z = P.Z * t0;           // ZP * (K1 * (XP + ZP) - K2 * (XP - ZP))^2
}

template <const mpz_class& p>
void MontgomeryCurve<p>::iso3e(const FpElem<p>& e, const MontgomeryPoint<p>& G0, MontgomeryCurve<p>& curve, 
                        const MontgomeryPoint<p>& P, const MontgomeryPoint<p>& Q, const MontgomeryPoint<p>& R,
                              MontgomeryPoint<p>& P1,      MontgomeryPoint<p>& Q1,      MontgomeryPoint<p>& R1) const
{
    MontgomeryPoint<p> G(G0);
    MontgomeryPoint<p> T;
    MontgomeryPoint<p> rG;
    MontgomeryPoint<p> rP1;
    MontgomeryPoint<p> rQ1;
    MontgomeryPoint<p> rR1;
    Fp2Elem<p> K1;
    Fp2Elem<p> K2;
    P1 = P.get_normalized();
    Q1 = Q.get_normalized();
    R1 = R.get_normalized();
    curve = *this;
    for (auto i = e - FpElem<p>("1"); i >= FpElem<p>("0"); --i)
    {
        curve.xTrpl_e(i, G, T);              // Ti = [3^i]Gi
        curve.iso3_curve(T, curve, K1, K2);  // Ei+1 = Ei/<Ti>
        curve.iso3_eval(T, P1, K1, K2, rP1); // Pi+1= phi3_i(Pi)
        curve.iso3_eval(T, Q1, K1, K2, rQ1); // Qi+1= phi3_i(Qi)
        curve.iso3_eval(T, R1, K1, K2, rR1); // Ri+1= phi3_i(Ri)
        P1 = rP1;
        Q1 = rQ1;
        R1 = rR1;
        if(i == FpElem<p>("0"))
            return;
        curve.iso3_eval(T, G, K1, K2, rG);   // Gi+1= phi3_i(Gi)
        G = rG;
    }
}

template <const mpz_class& p>
Fp2Elem<p> MontgomeryCurve<p>::j_inv(void) const
{
    Fp2Elem<p> j(A * A);    // A^2
    Fp2Elem<p> t1(C * C);   // C^2
    Fp2Elem<p> t0(t1 + t1); // 2C^2
    t0 = j - t0;            // A^2 - 2C^2
    t0 = t0 - t1;           // A^2 - 3C^2
    j = t0 - t1;            // A^2 - 4C^2
    t1 = t1 * t1;           // C^4
    j = j * t1;             // (A^2 - 4C^2) * C^4
    t0 = t0 + t0;           // 2*(A^2 - 3C^2)
    t0 = t0 + t0;           // 4*(A^2 - 3C^2)
    t1 = t0 * t0;           // (4*(A^2 - 3C^2))^2
    t0 = t0 * t1;           // 4*(A^2 - 3C^2) * (4*(A^2 - 3C^2))^2 = 4^3 *(A^2 - 3C^2)^3
    t0 = t0 + t0;           // 2 * (4^3 *(A^2 - 3C^2)^3)
    t0 = t0 + t0;           // 4 * (4^3 *(A^2 - 3C^2)^3)
    j = j.invert();         // 1 / ((A^2 - 4C^2) * C^4)
    j = j * t0;             // 4 * (4^3 *(A^2 - 3C^2)^3) / ((A^2 - 4C^2) * C^4)
    return j;               // 256 * (4^3 *(A^2 - 3C^2)^3) / ((A^2 - 4C^2) * C^4)
}

template <const mpz_class& p>
Fp2Elem<p> MontgomeryCurve<p>::get_A(void) const
{
    return A;
}

template <const mpz_class& p>
Fp2Elem<p> MontgomeryCurve<p>::get_C(void) const
{
    return C;
}

template <const mpz_class& p>
Fp2Elem<p> MontgomeryCurve<p>::get_a(void) const
{
    return A * C.invert();
}

template <const mpz_class& p>
MontgomeryCurve<p>::MontgomeryCurve(const MontgomeryPoint<p>& P, const MontgomeryPoint<p>& Q, const MontgomeryPoint<p>& R) : B(Fp2Elem<p>("1", "0")), C(Fp2Elem<p>("1", "0"))
{
    Fp2Elem<p> xP(P.get_x());
    Fp2Elem<p> xQ(Q.get_x());
    Fp2Elem<p> xR(R.get_x());

    Fp2Elem<p> t0(xP * xQ);       // xP * xQ
    Fp2Elem<p> t1(xP + xQ);       // xP + xQ
    A = xR * t1;                  // xR * (xP + xQ)
    A = A + t0;                   // xR * (xP + xQ) + xP * xQ
    t0 = t0 * xR;                 // xP * xQ * xR
    A = A - Fp2Elem<p>("1", "0"); // xPxR + xRxQ + xPxQ - 1
    t0 = t0 + t0;                 // 2 * xPxQxR
    t0 = t0 + t0;                 // 4 * xPxQxR
    t1 = t1 + xR;                 // xP + xQ + xR
    A = A * A;                    // (xPxR + xRxQ + xPxQ - 1)^2
    t0 = t0.invert();             // 1 / 4 * xPxQxR
    A = A * t0;                   // (xPxR + xRxQ + xPxQ - 1)^2 / 4xPxQxR
    A = A - t1;                   // (xPxR + xRxQ + xPxQ - 1)^2 / 4xPxQxR - (xP + xQ + xR)

    C24 = C + C;                  // 2C
    Ap24 = A + C24;               // A + 2C
    Am24 = A - C24;               // A - 2C
    C24 = C24 + C24;              // 4C
    // a24 = Ap24 * C.invert();      // A/C + 2
    Fp2Elem<p> t("2", "0");
    a24 = A + t;
    t = (t + t).invert();
    a24 = a24 * t;
}

template <const mpz_class& p>
void MontgomeryCurve<p>::print(void) const
{
    Fp2Elem<p> a = get_a();
    std::cout << "y^2 = x^3 + " << a.get_str()  << " * x^2 + x" << std::endl;
}