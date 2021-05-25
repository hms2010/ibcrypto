#pragma once
#include <gmpxx.h>
#include <iostream>
#include <string>

typedef mpz_class big_int_t;

template <const big_int_t& p>
class FpElem
{
        big_int_t _elem_;
    public:
        FpElem(void);
        ~FpElem(void) = default;

        FpElem(const big_int_t& elem);
        FpElem(const std::string& elem);

        FpElem(const FpElem& src);
        FpElem(FpElem&& src);


        FpElem& operator = (const FpElem& src);
        FpElem& operator = (FpElem&& src);

        FpElem operator + (void) const;
        FpElem operator + (const FpElem& op) const;
        FpElem operator + (FpElem&& op) const;
        FpElem& operator += (const FpElem& op);
        FpElem& operator += (FpElem&& op);

        FpElem operator - (void) const;
        FpElem operator - (const FpElem& op) const;
        FpElem operator - (FpElem&& op) const;
        FpElem& operator -= (const FpElem& op);
        FpElem& operator -= (FpElem&& op);

        FpElem operator * (const FpElem& op) const;
        FpElem operator * (FpElem&& op) const;
        FpElem& operator *= (const FpElem& op);
        FpElem& operator *= (FpElem&& op);

        FpElem operator / (const FpElem& op) const;
        FpElem operator / (FpElem&& op) const;
        FpElem& operator /= (const FpElem& op);
        FpElem& operator /= (FpElem&& op);

        FpElem operator % (const FpElem& op) const;
        FpElem operator % (FpElem&& op) const;
        FpElem& operator %= (const FpElem& op);
        FpElem& operator %= (FpElem&& op);

        FpElem invert(void);

        FpElem& operator ++ (void);
        FpElem& operator -- (void);

        bool operator == (const FpElem& op) const;
        bool operator != (const FpElem& op) const;
        bool operator < (const FpElem& op) const;
        bool operator <= (const FpElem& op) const;
        bool operator > (const FpElem& op) const;
        bool operator >= (const FpElem& op) const;

        std::string get_str(void) const;
        std::string get_p_str(void) const;
        void print(void) const;
};

template <const big_int_t& p>
class Fp2Elem
{
        FpElem<p> _real_;
        FpElem<p> _imag_;

    public:
        Fp2Elem(void);
        ~Fp2Elem(void) = default;
        Fp2Elem(const FpElem<p>& real, const FpElem<p>& imag);
        Fp2Elem(const std::string& real, const std::string& imag);
        Fp2Elem(const big_int_t& real, const big_int_t& imag);
        Fp2Elem(const Fp2Elem& src);
        Fp2Elem(Fp2Elem&& src) noexcept;

        Fp2Elem& operator = (const Fp2Elem& src);
        Fp2Elem& operator = (Fp2Elem&& src) noexcept;

        Fp2Elem operator + (void) const;
        Fp2Elem operator + (const Fp2Elem& op) const;
        Fp2Elem operator + (Fp2Elem&& op) const;
        Fp2Elem operator + (const FpElem<p>& op) const;
        Fp2Elem& operator += (const Fp2Elem& op);
        Fp2Elem& operator += (Fp2Elem&& op);

        Fp2Elem operator - (void) const;
        Fp2Elem operator - (const Fp2Elem& op) const;
        Fp2Elem operator - (Fp2Elem&& op) const;
        Fp2Elem& operator -= (const Fp2Elem& op);
        Fp2Elem& operator -= (Fp2Elem&& op);

        Fp2Elem operator * (const Fp2Elem& op) const;
        Fp2Elem operator * (Fp2Elem&& op) const;
        Fp2Elem& operator *= (const Fp2Elem& op);
        Fp2Elem& operator *= (Fp2Elem&& op);
        Fp2Elem invert(void) const;

        bool operator == (Fp2Elem& op);

        std::string get_real_str(void) const;
        std::string get_imag_str(void) const;
        std::string get_p_str(void) const;
        std::string get_str(void) const;
        void print(void) const;
};

#include "fp2_arithm_impl.cpp"
