#include <iostream>
#include <string>

#include "fp2_arithm.hpp"
#include "ec_arithm.hpp"

const mpz_class p = mpz_class("11");

bool test_FpElem(void)
{
    //static 
    FpElem<p> a("2");
    std::cout << "a = "; a.print();
    FpElem<p> b("9");
    std::cout << "b = "; b.print();
    FpElem<p> c = a + b;
    std::cout << "c = "; c.print();
    FpElem<p> d("19");
    std::cout << "d = "; d.print();
    c = d;
    std::cout << "c = "; c.print();
    d = FpElem<p>("8") - FpElem<p>("1");
    std::cout << "d = "; d.print();
    d = -d;
    std::cout << "d = "; d.print();
    a =  d * b;
    std::cout << "a = "; a.print();
    b = a.invert();
    std::cout << "b = "; b.print();
    FpElem<p> k = (a + b).invert();
    std::cout << "k = " << k.get_str() << " % " << k.get_p_str() << std::endl;
    k = k.invert();
    std::cout << "k = " << k.get_str() << " % " << k.get_p_str() << std::endl;
    std::cout << (a == a) << std::endl;
    std::cout << (a == k) << std::endl;
    return true;
}

bool test_Fp2Elem(void)
{
    Fp2Elem<p> a;
    std::cout << "a = "; a.print();
    Fp2Elem<p> b(mpz_class(1), mpz_class(4));
    std::cout << "b = "; b.print();
    Fp2Elem<p> c("3", "12");
    std::cout << "c = "; c.print();
    Fp2Elem<p> d = b + c;
    std::cout << "d = "; d.print();
    d = b - c;
    std::cout << "d = "; d.print();
    d = c - b;
    std::cout << "d = "; d.print();
    d = -d;
    std::cout << "d = "; d.print();
    c = d.invert();
    std::cout << "c = (" << d.get_real_str() << " + " << c.get_imag_str() << " * i) % " << d.get_p_str() << std::endl;
    (c * d).print();
    std::cout << (a == a) << std::endl;
    std::cout << (a == c) << std::endl;
    return true;
}

bool test_ec_arithm(void)
{
    Fp2Elem<p> x("9", "0");
    Fp2Elem<p> y("2", "0");
    MontgomeryCurve<p> mont_curve(x, y, Fp2Elem<p>("1", "0"));
    MontgomeryPoint<p> P(x, y);
    MontgomeryPoint<p> S(x, y);
    x = x + x;
    MontgomeryPoint<p> Q(x, y);
    FpElem<p> k("4");
    mont_curve.ladder3pt(k, P, Q, Q, S);

    MontgomeryPoint<p> T;
    mont_curve.xDblAdd(P, Q, Q, S, T);
    mont_curve.j_inv();

    MontgomeryPoint<p> O(Fp2Elem<p>("1", "0"), Fp2Elem<p>("0", "0"));
    mont_curve.xDblAdd(P, O, P, Q, S);
    O.print();
    P.print();
    P.get_normalized().print();
    Q.print();
    S.print();
    MontgomeryPoint<p> H(std::move(MontgomeryPoint<p>(P.X, P.Z)));
    std::cout << "H:";
    H.print();
    mont_curve.xDbl(P, Q);
    mont_curve.xDbl_e(k, P, Q);
    mont_curve.xTrpl(P, Q);
    MontgomeryCurve<p> iso_mont_curve(x, y, Fp2Elem<p>("1", "0"));
    mont_curve.iso2e(k, O, iso_mont_curve, P, Q, S, P, Q, S);
    MontgomeryCurve<p> from_point_mont_curve(P, Q, S);
    return true;
}

const mpz_class q = mpz_class("431");
bool test_proto(void)
{
    MontgomeryCurve<q> start_curve(Fp2Elem<q>("104", "0"), Fp2Elem<q>("1", "0"), Fp2Elem<q>("1", "0"));

    MontgomeryPoint<q> P2(Fp2Elem<q>("349", "381"));
    MontgomeryPoint<q> Q2(Fp2Elem<q>("52", "255"));
    MontgomeryPoint<q> R2(Fp2Elem<q>("357", "351"));

    MontgomeryPoint<q> P3(Fp2Elem<q>("126", "111"));
    MontgomeryPoint<q> Q3(Fp2Elem<q>("219", "199"));
    MontgomeryPoint<q> R3(Fp2Elem<q>("98", "14"));

    std::cout << "A = ";
    start_curve.get_A().print();
    std::cout << "C = ";
    start_curve.get_C().print();
    std::cout << "j_inv = ";
    start_curve.j_inv().print();

    FpElem<q> n2("11");
    FpElem<q> e2("4");
    MontgomeryCurve<q> curve2;
    MontgomeryCurve<q> curve32;
    MontgomeryPoint<q> G2;
    MontgomeryPoint<q> G32;
    MontgomeryPoint<q> P32;
    MontgomeryPoint<q> Q32;
    MontgomeryPoint<q> R32;

    FpElem<q> n3("19");
    FpElem<q> e3("3");
    MontgomeryCurve<q> curve3;
    MontgomeryCurve<q> curve23;
    MontgomeryPoint<q> G3;
    MontgomeryPoint<q> G23;
    MontgomeryPoint<q> P23;
    MontgomeryPoint<q> Q23;
    MontgomeryPoint<q> R23;


    start_curve.ladder3pt(n2, P2, Q2, R2, G2);
    std::cout << "G2 = ";
    G2.get_x().print();
    start_curve.iso2e(e2, G2, curve2, P3, Q3, R3, P23, Q23, R23);
    curve2.print();
    P23.get_x().print();
    Q23.get_x().print();
    R23.get_x().print();

    start_curve.ladder3pt(n3, P3, Q3, R3, G3);
    std::cout << "G3 = ";
    G3.get_x().print();
    start_curve.iso3e(e3, G3, curve3, P2, Q2, R2, P32, Q32, R32);
    curve3.print();
    P32.get_x().print();
    Q32.get_x().print();
    R32.get_x().print();

    MontgomeryPoint<q> P;
    curve3 = MontgomeryCurve<q>(P32, Q32, R32);
    curve3.print();
    curve3.ladder3pt(n2, P32, Q32, R32, G32);
    std::cout << "G32 = ";
    G32.get_x().print();
    curve3.iso2e(e2, G32, curve32, P32, Q32, R32, P, P, P);
    curve32.print();

    curve2 = MontgomeryCurve<q>(P23, Q23, R23);
    curve2.print();
    curve2.ladder3pt(n3, P23, Q23, R23, G23);
    std::cout << "G23 = ";
    G23.get_x().print();
    curve2.iso3e(e3, G23, curve23, P23, Q23, R23, P, P, P);
    curve23.print();

    Fp2Elem<q> j32(curve32.j_inv());
    Fp2Elem<q> j23(curve23.j_inv());

    bool res = j32 == j23;
    if(res)
        std::cout << "OK. j_inv =" << j23.get_str() << std::endl;
    else
        std::cout << "FAILED. j23_inv = " << j23.get_str() << ", j32_inv = " << j32.get_str() << std::endl;
    return res;
}

int main(void)
{
    // Fp2Elem<p> a("2", "1");
    // Fp2Elem<p> b(a.invert());
    // a.print();
    // b.print();
    // (a * b).print();
    // b = a * a;
    // b.print();
    // a *= a;
    // a.print();
    // test_FpElem();
    // test_Fp2Elem();
    // test_ec_arithm();
    test_proto();
    return 0;
}
