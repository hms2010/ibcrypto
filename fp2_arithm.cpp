#include "fp2_ariphm.hpp"

#include <iostream>

#include <gmp.h>
#include <gmpxx.h>

template <const mpz_class& p>
fp_elem<p>::fp_elem(const mpz_class& elem): _elem_(elem)
{
    _elem_ %= p;
}

template <const mpz_class& p>
fp_elem<p>::fp_elem(const std::string& elem): _elem_(elem)
{
    _elem_ %= p;
}

template <const mpz_class& p>
fp_elem<p>::fp_elem(const fp_elem<p>& src) :  _elem_(src._elem_)
{}

template <const mpz_class& p>
fp_elem<p>::fp_elem(fp_elem<p>&& src) : _elem_(std::move(src._elem_))
{}

template <const mpz_class& p>
fp_elem<p>& fp_elem<p>::operator = (const fp_elem<p>& src)
{
    _elem_ = src._elem_;
    return *this;
}

template <const mpz_class& p>
fp_elem<p>& fp_elem<p>::operator = (fp_elem<p>&& src)
{
    _elem_ = std::move(src._elem_);
    return *this;
}

template <const mpz_class& p>
fp_elem<p> fp_elem<p>::operator + (void) const
{
    return *this;
}

template <const mpz_class& p>
fp_elem<p> fp_elem<p>::operator + (const fp_elem<p>& op) const
{
    return fp_elem(_elem_ + op._elem_ % p);
}

template <const mpz_class& p>
fp_elem<p> fp_elem<p>::operator + (fp_elem<p>&& op) const
{
    return fp_elem(_elem_ + op._elem_ % p);
}

template <const mpz_class& p>
fp_elem<p>& fp_elem<p>::operator += (const fp_elem<p>& op)
{
    _elem_ = _elem_ + op._elem_ % p;
    return *this;
}

template <const mpz_class& p>
fp_elem<p>& fp_elem<p>::operator += (fp_elem<p>&& op)
{
    _elem_ = _elem_ + op._elem_ % p;
    return *this;
}

template <const mpz_class& p>
fp_elem<p> fp_elem<p>::operator - (void) const
{
    return fp_elem(p - _elem_);
}

template <const mpz_class& p>
fp_elem<p> fp_elem<p>::operator - (const fp_elem<p>& op) const
{
    return fp_elem(_elem_ + (p - op._elem_) % p);
}

template <const mpz_class& p>
fp_elem<p> fp_elem<p>::operator - (fp_elem<p>&& op) const
{
    return fp_elem(_elem_ + (p - op._elem_) % p);
}

template <const mpz_class& p>
fp_elem<p>& fp_elem<p>::operator -= (const fp_elem<p>& op)
{
    _elem_ = _elem_ + (p - op._elem_) % p;
    return *this;
}

template <const mpz_class& p>
fp_elem<p>& fp_elem<p>::operator -= (fp_elem<p>&& op)
{
    _elem_ = _elem_ + (p - op._elem_) % p;
    return *this;
}

template <const mpz_class& p>
fp_elem<p> fp_elem<p>::operator * (const fp_elem<p>& op) const
{
    return fp_elem(_elem_ * op._elem_ % p);
}

template <const mpz_class& p>
fp_elem<p> fp_elem<p>::operator * (fp_elem<p>&& op) const
{
    return fp_elem(_elem_ * op._elem_ % p);
}

template <const mpz_class& p>
fp_elem<p>& fp_elem<p>::operator *= (const fp_elem<p>& op)
{
    _elem_ = _elem_ * op._elem_ % p;
    return *this;
}

template <const mpz_class& p>
fp_elem<p>& fp_elem<p>::operator *= (fp_elem<p>&& op)
{
    _elem_ = _elem_ * op._elem_ % p;
    return *this;
}

template <const mpz_class& p>
fp_elem<p> fp_elem<p>::invert(void)
{
    mpz_t tmp;
    mpz_init(tmp);
    mpz_invert(tmp, _elem_.get_mpz_t(), p.get_mpz_t());
    fp_elem<p> res(mpz_class(tmp).get_str());
    mpz_clear(tmp);
    return res;
}

template <const mpz_class& p>
void fp_elem<p>::print(void) const
{
    std::cout << _elem_.get_str() << " % " << p.get_str() << std::endl;
}

template <const mpz_class& p>
std::string fp_elem<p>::get_str(void) const
{
    return _elem_.get_str();
}

template <const mpz_class& p>
std::string fp_elem<p>::get_p_str(void) const
{
    return p.get_str();
}

template <const mpz_class& p>
fp2_elem<p>::fp2_elem(const fp_elem<p>& real, const fp_elem<p>& imag) : _real_(real), _imag_(imag)
{}

// template <const mpz_class& p>
// fp2_elem<p>::fp2_elem (fp_elem<p>&& real, fp_elem<p>&& imag) : _real_(real), _imag_(imag)
// {}

template <const mpz_class& p>
fp2_elem<p>::fp2_elem(const mpz_class& real, const mpz_class& imag) : _real_(real), _imag_(imag)
{}

template <const mpz_class& p>
fp2_elem<p>::fp2_elem(const std::string& real, const std::string& imag) : _real_(real), _imag_(imag)
{}

template <const mpz_class& p>
fp2_elem<p>::fp2_elem(const fp2_elem<p>& src) : _real_(src._real_), _imag_(src._imag_)
{}

template <const mpz_class& p>
fp2_elem<p>::fp2_elem(fp2_elem<p>&& src) noexcept
{
    _real_ = std::move(src._real_);
    _imag_ = std::move(src._imag_);
}

template <const mpz_class& p>
fp2_elem<p>& fp2_elem<p>::operator = (const fp2_elem<p>& src)
{
    _real_ = src._real_;
    _imag_ = src._imag_;
    return *this;
}

template <const mpz_class& p>
fp2_elem<p>& fp2_elem<p>::operator = (fp2_elem<p>&& src) noexcept
{
    _real_ = std::move(src._real_);
    _imag_ = std::move(src._imag_);
    return *this;
}

template <const mpz_class& p>
fp2_elem<p> fp2_elem<p>::operator + (void) const
{
    return *this;
}

template <const mpz_class& p>
fp2_elem<p> fp2_elem<p>::operator + (const fp2_elem<p>& op) const
{
    return fp2_elem(_real_ + op._real_, _imag_ + op._imag_);
}

template <const mpz_class& p>
fp2_elem<p> fp2_elem<p>::operator + (fp2_elem<p>&& op) const
{
    return fp2_elem(_real_ + op._real_, _imag_ + op._imag_);
}

template <const mpz_class& p>
fp2_elem<p>& fp2_elem<p>::operator += (const fp2_elem<p>& op)
{
    _real_ = _real_ + op._real_;
    _imag_ = _imag_ + op._imag_;
    return *this;
}

template <const mpz_class& p>
fp2_elem<p>& fp2_elem<p>::operator += (fp2_elem<p>&& op)
{
    _real_ = _real_ + op._real_;
    _imag_ = _imag_ + op._imag_;
    return *this;
}

template <const mpz_class& p>
fp2_elem<p> fp2_elem<p>::operator - (void) const
{
    return fp2_elem(-_real_, _imag_);
}

template <const mpz_class& p>
fp2_elem<p> fp2_elem<p>::operator - (const fp2_elem<p>& op) const
{
    return fp2_elem(_real_ - op._real_, _imag_ - op._imag_);
}

template <const mpz_class& p>
fp2_elem<p> fp2_elem<p>::operator - (fp2_elem<p>&& op) const
{
    return fp2_elem(_real_ - op._real_, _imag_ - op._imag_);
}

template <const mpz_class& p>
fp2_elem<p>& fp2_elem<p>::operator -= (const fp2_elem<p>& op)
{
    _real_ = _real_ - op._real_;
    _imag_ = _imag_ - op._imag_;
    return *this;
}

template <const mpz_class& p>
fp2_elem<p>& fp2_elem<p>::operator -= (fp2_elem<p>&& op)
{
    _real_ = _real_ - op._real_;
    _imag_ = _imag_ - op._imag_;
    return *this;
}

template <const mpz_class& p>
fp2_elem<p> fp2_elem<p>::operator * (const fp2_elem<p>& op) const
{
    return fp2_elem(_real_ * op._real_ - _imag_ * op._imag_, _real_ * op._imag_ + _imag_ * op._real_);
}
    
template <const mpz_class& p>
fp2_elem<p> fp2_elem<p>::operator * (fp2_elem<p>&& op) const
{
    return fp2_elem(_real_ * op._real_ - _imag_ * op._imag_, _real_ * op._imag_ + _imag_ * op._real_);
}

template <const mpz_class& p>
fp2_elem<p>& fp2_elem<p>::operator *= (const fp2_elem<p>& op)
{
    _real_ = _real_ * op._real_ - _imag_ * op._imag_;
    _imag_ = _real_ * op._imag_ + _imag_ * op._real_;
    return *this;
}

template <const mpz_class& p>
fp2_elem<p>& fp2_elem<p>::operator *= (fp2_elem<p>&& op)
{
    _real_ = _real_ * op._real_ - _imag_ * op._imag_;
    _imag_ = _real_ * op._imag_ + _imag_ * op._real_;
    return *this;
}

template <const mpz_class& p>
fp2_elem<p> fp2_elem<p>::invert(void) const
{
    fp_elem<p> tmp = (_real_ * _real_ + _imag_ * _imag_).invert();
    fp2_elem<p> res;
    res._real_ = _real_ * tmp;
    res._imag_ = - _imag_;
    res._imag_ *= tmp;
    return res; 
}

template <const mpz_class& p>
std::string fp2_elem<p>::get_real_str(void) const
{
    return _real_.get_str();
}

template <const mpz_class& p>
std::string fp2_elem<p>::get_imag_str(void) const
{
    return _imag_.get_str();
}

template <const mpz_class& p>
std::string fp2_elem<p>::get_p_str(void) const
{
    return p.get_str();
}

template <const mpz_class& p>
void fp2_elem<p>::print(void) const
{
    std::cout << "(" <<_real_.get_str() << " + " << _imag_.get_str() << " * i) % " << p << std::endl;
}

#if defined(DEBUG)

const mpz_class p = mpz_class("11");
bool test_fp_elem(void)
{
    //static 
    fp_elem<p> a("2");
    std::cout << "a = "; a.print();
    fp_elem<p> b("9");
    std::cout << "b = "; b.print();
    fp_elem<p> c = a + b;
    std::cout << "c = "; c.print();
    fp_elem<p> d("19");
    std::cout << "d = "; d.print();
    c = d;
    std::cout << "c = "; c.print();
    d = fp_elem<p>("8") - fp_elem<p>("1");
    std::cout << "d = "; d.print();
    d = -d;
    std::cout << "d = "; d.print();
    a =  d * b;
    std::cout << "a = "; a.print();
    b = a.invert();
    std::cout << "b = "; b.print();
    fp_elem<p> k = (a + b).invert();
    std::cout << "k = " << k.get_str() << " % " << k.get_p_str() << std::endl;
    k = k.invert();
    std::cout << "k = " << k.get_str() << " % " << k.get_p_str() << std::endl;
    return true;
}

bool test_fp2_elem(void)
{
    fp2_elem<p> a;
    std::cout << "a = "; a.print();
    fp2_elem<p> b(mpz_class(1), mpz_class(4));
    std::cout << "b = "; b.print();
    fp2_elem<p> c("3", "12");
    std::cout << "c = "; c.print();
    fp2_elem<p> d = b + c;
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
    // fp2_elem<p> k = -b;
    // b = -b;
    return true;
}

int main(void)
{
    test_fp_elem();
    test_fp2_elem();
    return 0;
}
#endif // DEBUG
