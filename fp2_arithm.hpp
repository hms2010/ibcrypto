#pragma once
#include <gmpxx.h>

typedef mpz_class big_int_t;

// template <big_int_t& _p_>
class fp_elem
{
        big_int_t _elem_;
        // big_int_t _p_;
    public:
        fp_elem(void)=delete;
        ~fp_elem(void)=default;
        fp_elem(const big_int_t& p);
        fp_elem(const big_int_t& elem, const big_int_t& p);
        fp_elem(const std::string& p);
        fp_elem(const std::string& elem, const std::string& p);
        fp_elem(const fp_elem& src);
        fp_elem(fp_elem&& src);

        fp_elem& operator = (const fp_elem& src);
        fp_elem& operator = (fp_elem&& src);

        fp_elem operator + (void) const;
        fp_elem operator + (const fp_elem& op) const;
        fp_elem operator + (fp_elem&& op) const;
        fp_elem& operator += (const fp_elem& op);
        fp_elem& operator += (fp_elem&& op);

        fp_elem operator - (void);
        fp_elem operator - (const fp_elem& op) const;
        fp_elem operator - (fp_elem&& op) const;
        fp_elem& operator -= (const fp_elem&);
        fp_elem& operator -= (fp_elem&&);

        fp_elem operator * (const fp_elem& op) const;
        fp_elem operator * (fp_elem&& op) const;
        fp_elem& operator *= (const fp_elem&);
        fp_elem& operator *= (fp_elem&&);

        friend fp_elem invert(const fp_elem& op);
        friend fp_elem invert(fp_elem&& op);
};

class fp2_elem
{
        fp_elem _real_;
        fp_elem _imag_;

    public:
        friend class fp2_field;
        fp2_elem(void)=delete;
        ~fp2_elem(void)=default;
        fp2_elem(const fp2_elem&);
        fp2_elem(fp2_elem&&) noexcept;

        fp2_elem& operator = (const fp2_elem&);
        fp2_elem& operator = (fp2_elem&&) noexcept;

        fp2_elem operator + (void) const;
        fp2_elem operator + (const fp2_elem&) const;
        fp2_elem operator + (fp2_elem&&) const;
        fp2_elem operator + (const fp_elem&) const;
        friend fp2_elem operator + (const fp_elem&, const fp2_elem&);
        fp2_elem& operator += (const fp2_elem&);
        fp2_elem& operator += (fp2_elem&&);

        fp2_elem operator - (void) const;
        fp2_elem operator - (const fp2_elem&) const;
        fp2_elem operator - (fp2_elem&&) const;
        friend fp2_elem operator - (const fp_elem&, const fp2_elem&);
        fp2_elem& operator -= (const fp2_elem&);
        fp2_elem& operator -= (fp2_elem&&);

        fp2_elem operator * (const fp2_elem&) const;
        fp2_elem operator * (fp2_elem&&) const;
        friend fp2_elem operator * (const fp_elem&, const fp2_elem&);
        fp2_elem& operator *= (const fp2_elem&);
        fp2_elem& operator *= (fp2_elem&&);

        friend fp2_elem invert(const fp2_elem&);
        friend fp2_elem invert(fp2_elem&&);
        fp2_elem operator / (const fp2_elem&) const;
        fp2_elem operator / (fp2_elem&&) const;

        friend fp2_elem operator / (const fp_elem&, const fp2_elem&);
        fp2_elem& operator /= (const fp2_elem&);
        fp2_elem& operator /= (fp2_elem&&);
};

class fp2_field
{
        big_int_t p;
    public:
        fp2_field(big_int_t&);
        fp2_field(std::string);
        ~fp2_field(void)=default;

        fp2_field(void)=delete;
        fp2_field(const fp2_field&)=delete;
        fp2_field(fp2_field&&)=delete;
};
