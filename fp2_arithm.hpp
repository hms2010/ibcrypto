#pragma once
#include <gmpxx.h>
#include <iostream>

#define GMP_ARITHM_IMPL

#ifdef GMP_ARITHM_IMPL
typedef mpz_class big_int_t;
#else
#error F_p ariphmetic implemenation needs to be defined
#endif

template <const big_int_t& p>
class fp_elem
{
        big_int_t _elem_;
    public:
        fp_elem(void)=default;
        ~fp_elem(void)=default;
        fp_elem(const big_int_t& elem);
        fp_elem(const std::string& elem);
        fp_elem(const fp_elem& src);
        fp_elem(fp_elem&& src);

        fp_elem& operator = (const fp_elem& src);
        fp_elem& operator = (fp_elem&& src);

        fp_elem operator + (void) const;
        fp_elem operator + (const fp_elem& op) const;
        fp_elem operator + (fp_elem&& op) const;
        fp_elem& operator += (const fp_elem& op);
        fp_elem& operator += (fp_elem&& op);

        fp_elem operator - (void) const;
        fp_elem operator - (const fp_elem& op) const;
        fp_elem operator - (fp_elem&& op) const;
        fp_elem& operator -= (const fp_elem& op);
        fp_elem& operator -= (fp_elem&& op);

        fp_elem operator * (const fp_elem& op) const;
        fp_elem operator * (fp_elem&& op) const;
        fp_elem& operator *= (const fp_elem& op);
        fp_elem& operator *= (fp_elem&& op);
        fp_elem invert(void);

        std::string get_str(void) const;
        std::string get_p_str(void) const;
        void print(void) const;
};

template <const big_int_t& p>
class fp2_elem
{
        fp_elem<p> _real_;
        fp_elem<p> _imag_;

    public:
        fp2_elem(void)=default;
        ~fp2_elem(void)=default;
        fp2_elem(const fp_elem<p>& real, const fp_elem<p>& imag);
        fp2_elem(const std::string& real, const std::string& imag);
        fp2_elem(const big_int_t& real, const big_int_t& imag);
        fp2_elem(const fp2_elem& src);
        fp2_elem(fp2_elem&& src) noexcept;

        fp2_elem& operator = (const fp2_elem& src);
        fp2_elem& operator = (fp2_elem&& src) noexcept;

        fp2_elem operator + (void) const;
        fp2_elem operator + (const fp2_elem& op) const;
        fp2_elem operator + (fp2_elem&& op) const;
        fp2_elem operator + (const fp_elem<p>& op) const;
        fp2_elem& operator += (const fp2_elem& op);
        fp2_elem& operator += (fp2_elem&& op);

        fp2_elem operator - (void) const;
        fp2_elem operator - (const fp2_elem& op) const;
        fp2_elem operator - (fp2_elem&& op) const;
        fp2_elem& operator -= (const fp2_elem& op);
        fp2_elem& operator -= (fp2_elem&& op);

        fp2_elem operator * (const fp2_elem& op) const;
        fp2_elem operator * (fp2_elem&& op) const;
        fp2_elem& operator *= (const fp2_elem& op);
        fp2_elem& operator *= (fp2_elem&& op);

        fp2_elem invert(void) const;
        std::string get_real_str(void) const;
        std::string get_imag_str(void) const;
        std::string get_p_str(void) const;
        void print(void) const;
};
