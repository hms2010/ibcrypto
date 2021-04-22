#include "fp2_ariphm.hpp"

#if defined(DEBUG)
#include <iostream>
#endif // DEBUG

#include <gmp.h>
#include <gmpxx.h>
typedef mpz_class big_int_t;

// fp_elem::fp_elem(const big_int_t& p) : _elem_(0), _p_(p)
// {}

fp_elem::fp_elem(const big_int_t& elem, const big_int_t& p) : _elem_(elem), _p_(p)
{}

fp_elem::fp_elem(const std::string& p) : _elem_(0), _p_(p)
{}

fp_elem::fp_elem(const std::string& elem, const std::string& p): _elem_(elem), _p_(p)
{}

fp_elem::fp_elem(const fp_elem& src) :  _elem_(src._elem_), _p_(src._p_)
{}

fp_elem::fp_elem(fp_elem&& src) : _elem_(std::move(src._elem_)), _p_(std::move(src._p_))
{}

fp_elem& fp_elem::operator = (const fp_elem& src)
{
    _elem_ = src._elem_;
    _p_ = src._p_;
    return *this;
}

fp_elem& fp_elem::operator = (fp_elem&& src)
{
    _elem_ = std::move(src._elem_);
    _p_ = std::move(src._p_);
    return *this;
}

fp_elem fp_elem::operator + (void) const
{
    return *this;
}

fp_elem fp_elem::operator + (const fp_elem& op) const
{
    if(_p_ != op._p_)
        ; // throw exception here
    return fp_elem(_elem_ + op._elem_ % _p_);
}

fp_elem fp_elem::operator + (fp_elem&& op) const
{
    if(_p_ != op._p_)
        ; // throw exception here
    return fp_elem(_elem_ + op._elem_ % _p_);
}

fp_elem& fp_elem::operator += (const fp_elem& op)
{
    if(_p_ != op._p_)
        ; // throw exception here
    _elem_ = _elem_ + op._elem_ % _p_;
    return *this;
}

fp_elem& fp_elem::operator += (fp_elem&& op)
{
    if(_p_ != op._p_)
        ; // throw exception here
    _elem_ = _elem_ + op._elem_ % _p_;
    return *this;
}

fp_elem fp_elem::operator - (void)
{
    return fp_elem(_p_ - _elem_);
}

fp_elem fp_elem::operator - (const fp_elem& op) const
{
    if(_p_ != op._p_)
        ; // throw exception here
    return fp_elem(_elem_ + (_p_ - op._elem_) % _p_);
}

fp_elem fp_elem::operator - (fp_elem&& op) const
{
    if(_p_ != op._p_)
        ; // throw exception here
    return fp_elem(_elem_ + (_p_ - op._elem_) % _p_);
}

fp_elem& fp_elem::operator -= (const fp_elem& op)
{
    if(_p_ != op._p_)
        ; // throw exception here
    _elem_ = _elem_ + (_p_ - op._elem_) % _p_;
    return *this;
}

fp_elem& fp_elem::operator -= (fp_elem&& op)
{
    if(_p_ != op._p_)
        ; // throw exception here
    _elem_ = _elem_ + (_p_ - op._elem_) % _p_;
    return *this;
}

fp_elem fp_elem::operator * (const fp_elem& op) const
{
    if(_p_ != op._p_)
        ; // throw exception here
    return fp_elem(_elem_ * op._elem_ % _p_);
}

fp_elem fp_elem::operator * (fp_elem&& op) const
{
    if(_p_ != op._p_)
        ; // throw exception here
    return fp_elem(_elem_ * op._elem_ % _p_);
}

fp_elem& fp_elem::operator *= (const fp_elem& op)
{
    if(_p_ != op._p_)
        ; // throw exception here
    _elem_ = _elem_ * op._elem_ % _p_;
    return *this;
}

fp_elem& fp_elem::operator *= (fp_elem&& op)
{
    if(_p_ != op._p_)
        ; // throw exception here
    _elem_ = _elem_ * op._elem_ % _p_;
    return *this;
}

fp_elem invert(const fp_elem& op)
{
    mpz_t tmp;
    mpz_invert(tmp, op._elem_.get_mpz_t(), op._p_.get_mpz_t());
    fp_elem res(mpz_class(tmp), op._p_);
    mpz_clear(tmp);
    return res;
}

fp_elem invert(fp_elem&& op)
{
    mpz_t tmp;
    mpz_invert(tmp, op._elem_.get_mpz_t(), op._p_.get_mpz_t());
    fp_elem res(mpz_class(tmp), op._p_);
    mpz_clear(tmp);
    return res;
}


#if 0
fp2_elem::fp2_elem(mpz_class& _p) : p(_p), _real_(0), _imag_(0)
{
    #if defined(DEBUG)
    std::cout << "ctor default: " << _real_ << " + " << _imag_ << " * i" << std::endl; 
    #endif // DEBUG
}

fp2_elem::fp2_elem(const mpz_class& _p, const mpz_class& __real_, const mpz_class& __imag_) : p(_p), _real_(__real_ % p), _imag_(__imag_ % p)
{
    #if defined(DEBUG)
    std::cout << "ctor: " << _real_ << " + " << _imag_ << " * i" << std::endl; 
    #endif // DEBUG
}

fp2_elem::fp2_elem(const fp2_elem& src) : _real_(src._real_ % p), _imag_(src._imag_ % p)
{
    #if defined(DEBUG)
    std::cout << "ctor copy: " << _real_ << " + " << _imag_ << " * i" << std::endl; 
    #endif // DEBUG
}
fp2_elem::fp2_elem(fp2_elem&& src) noexcept
{
    _real_ = std::move(src._real_);
    _imag_ = std::move(src._imag_);
    #if defined(DEBUG)
    std::cout << "ctor move: " << _real_ << " + " << _imag_ << " * i" << std::endl; 
    #endif // DEBUG
}

fp2_elem& fp2_elem::operator = (const fp2_elem& src)
{
    #if defined(DEBUG)
    #endif // DEBUG
    _real_ = src._real_;
    _imag_ = src._imag_;
    std::cout << "fp2_elem& operator =: " << _real_ << " + " << _imag_ << " * i" << std::endl; 
    return *this;
}

fp2_elem& fp2_elem::operator = (fp2_elem&& src) noexcept
{
    #if defined(DEBUG)
    #endif // DEBUG
    _real_ = std::move(src._real_ % p);
    _imag_ = std::move(src._imag_ % p);
    std::cout << "fp2_elem&& operator =: " << _real_ << " + " << _imag_ << " * i" << std::endl; 
    return *this;
}

fp2_elem fp2_elem::operator + (void) const
{
    return *this;
}

fp2_elem fp2_elem::operator + (const fp2_elem& op) const
{
    return fp2_elem((_real_ + op._real_) % p, (_imag_ + op._imag_) % p);
}

fp2_elem fp2_elem::operator + (fp2_elem&& op) const
{
    return fp2_elem((_real_ + op._real_) % p, (_imag_ + op._imag_) % p);
}

fp2_elem fp2_elem::operator + (const mpz_class& op) const
{
    return fp2_elem((_real_ + op) % p, _imag_);
}
fp2_elem& fp2_elem::operator += (const fp2_elem& op)
{
    _real_ = (_real_ + op._real_) % p;
    _imag_ = (_imag_ + op._imag_) % p;
    return *this;
}

fp2_elem& fp2_elem::operator += (fp2_elem&& op)
{
    _real_ = (_real_ + op._real_) % p;
    _imag_ = (_imag_ + op._imag_) % p;
    return *this;
}


fp2_elem& fp2_elem::operator += (const mpz_class& op)
{
    _real_ = (_real_ + op) % p;
    return *this;
}

fp2_elem fp2_elem::operator - (void) const
{
    return fp2_elem(p - _real_ % p, p - _imag_ % p);
}

fp2_elem fp2_elem::operator - (const fp2_elem& op) const
{
    return fp2_elem((_real_ + (p - op._real_ % p)) % p, (_imag_ + (p - op._imag_ % p)) % p);
}

fp2_elem fp2_elem::operator - (fp2_elem&& op) const
{
    return fp2_elem((_real_ + (p - op._real_ % p)) % p, (_imag_ + (p - op._imag_ % p)) % p);
}

fp2_elem fp2_elem::operator - (const mpz_class& op) const
{
    return fp2_elem((_real_ + (p - op % p)) % p, _imag_ % p);
}
fp2_elem& fp2_elem::operator -= (const fp2_elem& op)
{
    _real_ = (_real_ + (p - op._real_ % p)) % p;
    _imag_ = (_imag_ + (p - op._imag_ % p)) % p;
    return *this;
}

fp2_elem& fp2_elem::operator -= (fp2_elem&& op)
{
    _real_ = (_real_ + (p - op._real_ % p)) % p;
    _imag_ = (_imag_ + (p - op._imag_ % p)) % p;
    return *this;
}


fp2_elem& fp2_elem::operator -= (const mpz_class& op)
{
    _real_ = (_real_ + (p - op % p)) % p;
    return *this;
}

fp2_elem fp2_elem::operator * (const fp2_elem& op) const
{
    return fp2_elem((_real_ * op._real_ % p + (p - _imag_ * op._imag_ % p)) % p, (_real_ * op._imag_ % p + _imag_ * op._real_ % p) % p);
}
    
fp2_elem fp2_elem::operator * (fp2_elem&& op) const
{
    return fp2_elem((_real_ * op._real_ % p + (p - _imag_ * op._imag_ % p)) % p, (_real_ * op._imag_ % p + _imag_ * op._real_ % p) % p);
}

fp2_elem fp2_elem::operator * (const mpz_class& op) const
{
    return fp2_elem(_real_ * op % p, _imag_ * op % p);
}

fp2_elem& fp2_elem::operator *= (const fp2_elem& op)
{
    _real_ = (_real_ * op._real_ % p + (p - _imag_ * op._imag_ % p)) % p;
    _imag_ = (_real_ * op._imag_ % p + _imag_ * op._real_) % p;
    return *this;
}

fp2_elem& fp2_elem::operator *= (fp2_elem&& op)
{
    _real_ = (_real_ * op._real_ % p + (p - _imag_ * op._imag_ % p)) % p;
    _imag_ = (_real_ * op._imag_ % p + _imag_ * op._real_) % p;
    return *this;
}

fp2_elem invert(const fp2_elem& op)
{
}

fp2_elem invert(fp2_elem&& op)
{
}

fp2_elem fp2_elem::operator / (const fp2_elem& op) const
{
    mpz_class denominator = op._real_ * op._real_ + op._imag_ * op._imag_;
    return fp2_elem((_real_ * op._real_ % p + (p - _imag_ * op._imag_ % p)) % p,
                    (_real_ * op._imag_ % p + _imag_ * op._real_ % p) % p);
}
    
fp2_elem fp2_elem::operator / (fp2_elem&& op) const
{
    mpz_class denominator = op._real_ * op._real_ + op._imag_ * op._imag_;
    return fp2_elem((_real_ * op._real_ + _imag_ * op._imag_) / denominator,
                 (-_real_ * op._imag_ + _imag_ * op._real_) / denominator);
}

fp2_elem fp2_elem::operator / (const mpz_class& op) const
{
    return fp2_elem(_real_ / op, _imag_ / op);
}

fp2_elem& fp2_elem::operator /= (const fp2_elem& op)
{
    mpz_class denominator = op._real_ * op._real_ + op._imag_ * op._imag_;
    _real_ = (_real_ * op._real_ + _imag_ * op._imag_) / denominator;
    _imag_ = (-_real_ * op._imag_ + _imag_ * op._real_) / denominator;
    return *this;
}

fp2_elem& fp2_elem::operator /= (fp2_elem&& op)
{
    mpz_class denominator = op._real_ * op._real_ + op._imag_ * op._imag_;
    _real_ = (_real_ * op._real_ + _imag_ * op._imag_) / denominator;
    _imag_ = (-_real_ * op._imag_ + _imag_ * op._real_) / denominator;
    return *this;
}

fp2_elem& fp2_elem:: operator /= (const mpz_class& op)
{
    _real_ = _real_ / op;
    _imag_ = _imag_ / op;
    return *this;
}
namespace fp2_field
{
    fp2_field(mpz_class& _p) : p(_p)
    {}

    fp2_elem& init_elem(void) const
    {
        return fp2_elem(p);
    }

    fp2_elem& init_elem(mpz_class& _real_, mpz_class& _imag_) const
    {
        return fp2_elem(p, _real_, _imag_);
    }
}
#endif // 0

#if defined(DEBUG)
bool test_fp_elem(void)
{
    fp_elem a()
}

int main(void)
{
#if 0
    fp2_elem a;
    fp2_elem b(mpz_class(1), mpz_class(4));
    fp2_elem c(mpz_class(3), mpz_class(2));
    fp2_elem d = b + c;
    d = d / 2;
    a = +b;
    b = -b;
#endif // 0

    return 0;
}
#endif // DEBUG
