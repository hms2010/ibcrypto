#include <iostream>
#include <string>
#include <ctime>
#include <string>

#include <gmp.h>
#include "fp2_arithm.hpp"
#include "ec_arithm.hpp"
#include "forsythia.hpp"

bool test_proto_toy(void);
bool test_proto(ForsythiaSet param_set, unsigned long& delta);

mpz_class p;
int main(void)
{
    long double common_time = 0;
    unsigned long delta = 0;
    long long runs = 100;
    bool res;

    p = mpz_class(g_param_set[forsythia128].p);
    for (long long e = 0; e < runs; ++e)
    {
        res = test_proto(forsythia128, delta);
        if(!res)
        {
            std::cout << "FAILED" << std::endl;
            break;
        }
        common_time += delta;
    }
    std::cout << "forsythia128: " << common_time / runs << "s" << std::endl;

    common_time = 0;
    p = mpz_class(g_param_set[forsythia256].p);
    for (long long e = 0; e < runs; ++e)
    {
        res = test_proto(forsythia256, delta);
        if(!res)
        {
            std::cout << "FAILED" << std::endl;
            break;
        }
        common_time += delta;
    }
    std::cout << "forsythia256: " << common_time / runs << "s" << std::endl;
    return 0;
}

bool test_proto(ForsythiaSet param_set, unsigned long & delta)
{
    std::time_t start_time;
    std::time_t end_time;

    Forsythia<p> for2(param_set, for_side_2);
    FpElem<p> sk2;
    MontgomeryPoint<p> P23;
    MontgomeryPoint<p> Q23;
    MontgomeryPoint<p> R23;
    Fp2Elem<p> j32;

    Forsythia<p> for3(param_set, for_side_3);
    FpElem<p> sk3;
    MontgomeryPoint<p> P32;
    MontgomeryPoint<p> Q32;
    MontgomeryPoint<p> R32;
    Fp2Elem<p> j23;
    start_time = std::time(nullptr);

    sk2.set_random(for2.get_e());
    for2.isogen(sk2, P23, Q23, R23);

    sk3.set_random(for3.get_e());
    for3.isogen(sk3, P32, Q32, R32);

    for2.isoex(sk2, P32, Q32, R32, j32);
    for3.isoex(sk3, P23, Q23, R23, j23);

    end_time = std::time(nullptr);
    delta = end_time - start_time;

    Fp2Elem<p> null_elem("0", "0");
    bool res = (j32 == j23) && !(j23 == null_elem);
    if(!res)
        std::cout << "FAILED. j23 = " << j23.get_str() << ", j32 = " << j32.get_str() << std::endl;
    #ifdef DEBUG
    else
        std::cout << "OK. j_inv = " << std::endl << j23.get_str() << std::endl;
    #endif // DEBUG
    return res;
}

const mpz_class q = mpz_class("431");
bool test_proto_toy(void)
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
