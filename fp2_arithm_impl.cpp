#include <iostream>
#include <ctime>

#include <gmp.h>
#include <gmpxx.h>

template <const mpz_class& p>
FpElem<p>::FpElem(void): _elem_("0")
{}

template <const mpz_class& p>
FpElem<p>::FpElem(const mpz_class& elem): _elem_(elem)
{
    _elem_ %= p;
}

template <const mpz_class& p>
FpElem<p>::FpElem(const std::string& elem): _elem_(elem)
{
    _elem_ %= p;
}

template <const mpz_class& p>
FpElem<p>::FpElem(const FpElem<p>& src) :  _elem_(src._elem_)
{}

template <const mpz_class& p>
FpElem<p>::FpElem(FpElem<p>&& src) : _elem_(std::move(src._elem_))
{}

template <const mpz_class& p>
FpElem<p>& FpElem<p>::operator = (const FpElem<p>& src)
{
    _elem_ = src._elem_;
    return *this;
}

template <const mpz_class& p>
FpElem<p>& FpElem<p>::operator = (FpElem<p>&& src)
{
    _elem_ = std::move(src._elem_);
    return *this;
}

template <const mpz_class& p>
FpElem<p> FpElem<p>::operator + (void) const
{
    return *this;
}

template <const mpz_class& p>
FpElem<p> FpElem<p>::operator + (const FpElem<p>& op) const
{
    return FpElem(_elem_ + op._elem_ % p);
}

template <const mpz_class& p>
FpElem<p> FpElem<p>::operator + (FpElem<p>&& op) const
{
    return FpElem(_elem_ + op._elem_ % p);
}

template <const mpz_class& p>
FpElem<p>& FpElem<p>::operator += (const FpElem<p>& op)
{
    _elem_ = _elem_ + op._elem_ % p;
    return *this;
}

template <const mpz_class& p>
FpElem<p>& FpElem<p>::operator += (FpElem<p>&& op)
{
    _elem_ = _elem_ + op._elem_ % p;
    return *this;
}

template <const mpz_class& p>
FpElem<p> FpElem<p>::operator - (void) const
{
    return FpElem(p - _elem_);
}

template <const mpz_class& p>
FpElem<p> FpElem<p>::operator - (const FpElem<p>& op) const
{
    return FpElem(_elem_ + (p - op._elem_) % p);
}

template <const mpz_class& p>
FpElem<p> FpElem<p>::operator - (FpElem<p>&& op) const
{
    return FpElem(_elem_ + (p - op._elem_) % p);
}

template <const mpz_class& p>
FpElem<p>& FpElem<p>::operator -= (const FpElem<p>& op)
{
    _elem_ = _elem_ + (p - op._elem_) % p;
    return *this;
}

template <const mpz_class& p>
FpElem<p>& FpElem<p>::operator -= (FpElem<p>&& op)
{
    _elem_ = _elem_ + (p - op._elem_) % p;
    return *this;
}

template <const mpz_class& p>
FpElem<p> FpElem<p>::operator * (const FpElem<p>& op) const
{
    return FpElem(_elem_ * op._elem_ % p);
}

template <const mpz_class& p>
FpElem<p> FpElem<p>::operator * (FpElem<p>&& op) const
{
    return FpElem(_elem_ * op._elem_ % p);
}

template <const mpz_class& p>
FpElem<p>& FpElem<p>::operator *= (const FpElem<p>& op)
{
    _elem_ = _elem_ * op._elem_ % p;
    return *this;
}

template <const mpz_class& p>
FpElem<p>& FpElem<p>::operator *= (FpElem<p>&& op)
{
    _elem_ = _elem_ * op._elem_ % p;
    return *this;
}

template <const mpz_class& p>
FpElem<p> FpElem<p>::operator / (const FpElem<p>& op) const
{
    return FpElem(_elem_ / op._elem_);
}

template <const mpz_class& p>
FpElem<p> FpElem<p>::operator / (FpElem<p>&& op) const
{
    return FpElem(_elem_ / op._elem_);
}

template <const mpz_class& p>
FpElem<p>& FpElem<p>::operator /= (const FpElem<p>& op)
{
    _elem_ = _elem_ / op._elem_;
    return *this;
}

template <const mpz_class& p>
FpElem<p>& FpElem<p>::operator /= (FpElem<p>&& op)
{
    _elem_ = _elem_ / op._elem_;
    return *this;
}

template <const mpz_class& p>
FpElem<p> FpElem<p>::operator % (const FpElem<p>& op) const
{
    return FpElem(_elem_ % op._elem_);
}

template <const mpz_class& p>
FpElem<p> FpElem<p>::operator % (FpElem<p>&& op) const
{
    return FpElem(_elem_ % op._elem_);
}

template <const mpz_class& p>
FpElem<p>& FpElem<p>::operator %= (const FpElem<p>& op)
{
    _elem_ = _elem_ % op._elem_;
    return *this;
}

template <const mpz_class& p>
FpElem<p>& FpElem<p>::operator %= (FpElem<p>&& op)
{
    _elem_ = _elem_ % op._elem_;
    return *this;
}

template <const mpz_class& p>
FpElem<p> FpElem<p>::invert(void)
{
    mpz_t tmp;
    mpz_init(tmp);
    mpz_invert(tmp, _elem_.get_mpz_t(), p.get_mpz_t());
    FpElem<p> res(mpz_class(tmp).get_str());
    mpz_clear(tmp);
    return res;
}

template <const mpz_class& p>
FpElem<p>& FpElem<p>::operator ++ (void)
{
    if (_elem_ == p - 1)
        _elem_ = 0; 
    else
        ++_elem_;
    return *this;
}

template <const mpz_class& p>
FpElem<p>& FpElem<p>::operator -- (void)
{
    if (_elem_ == 0)
        _elem_ = p - 1; 
    else
        --_elem_;
    return *this;
}

template <const mpz_class& p>
bool FpElem<p>::operator == (const FpElem<p>& op) const
{
    return cmp(_elem_, op._elem_) == 0;
}

template <const mpz_class& p>
bool FpElem<p>::operator != (const FpElem<p>& op) const
{
    return cmp(_elem_, op._elem_) != 0;
}

template <const mpz_class& p>
bool FpElem<p>::operator < (const FpElem<p>& op) const
{
    return cmp(_elem_, op._elem_) < 0;
}

template <const mpz_class& p>
bool FpElem<p>::operator <= (const FpElem<p>& op) const
{
    return cmp(_elem_, op._elem_) <= 0;
}

template <const mpz_class& p>
bool FpElem<p>::operator > (const FpElem<p>& op) const
{
    return cmp(_elem_, op._elem_) > 0;
}

template <const mpz_class& p>
bool FpElem<p>::operator >= (const FpElem<p>& op) const
{
    return cmp(_elem_, op._elem_) >= 0;
}

template <const mpz_class& p>
void FpElem<p>::print(void) const
{
    std::cout << _elem_.get_str() << " % " << p.get_str() << std::endl;
}

template <const mpz_class& p>
std::string FpElem<p>::get_str(void) const
{
    return _elem_.get_str();
}

template <const mpz_class& p>
std::string FpElem<p>::get_p_str(void) const
{
    return p.get_str();
}

template <const mpz_class& p>
void FpElem<p>::set_random(const unsigned long bitcnt)
{
    gmp_randclass random(gmp_randinit_default);
    random.seed(std::time(nullptr));
    _elem_ = random.get_z_bits(bitcnt);
}

template <const mpz_class& p>
unsigned long FpElem<p>::get_ui(void) const
{
    return mpz_get_ui(_elem_.get_mpz_t());
} 

template <const mpz_class& p>
Fp2Elem<p>::Fp2Elem(void) : _real_(FpElem<p>("0")), _imag_(FpElem<p>("0"))
{}

template <const mpz_class& p>
Fp2Elem<p>::Fp2Elem(const FpElem<p>& real, const FpElem<p>& imag) : _real_(real), _imag_(imag)
{}

template <const mpz_class& p>
Fp2Elem<p>::Fp2Elem(const mpz_class& real, const mpz_class& imag) : _real_(real), _imag_(imag)
{}

template <const mpz_class& p>
Fp2Elem<p>::Fp2Elem(const std::string& real, const std::string& imag) : _real_(real), _imag_(imag)
{}

template <const mpz_class& p>
Fp2Elem<p>::Fp2Elem(const Fp2Elem<p>& src) : _real_(src._real_), _imag_(src._imag_)
{}

template <const mpz_class& p>
Fp2Elem<p>::Fp2Elem(Fp2Elem<p>&& src) noexcept : _real_(std::move(src._real_)), _imag_(std::move(src._imag_))
{}

template <const mpz_class& p>
Fp2Elem<p>& Fp2Elem<p>::operator = (const Fp2Elem<p>& src)
{
    _real_ = src._real_;
    _imag_ = src._imag_;
    return *this;
}

template <const mpz_class& p>
Fp2Elem<p>& Fp2Elem<p>::operator = (Fp2Elem<p>&& src) noexcept
{
    _real_ = std::move(src._real_);
    _imag_ = std::move(src._imag_);
    return *this;
}

template <const mpz_class& p>
Fp2Elem<p> Fp2Elem<p>::operator + (void) const
{
    return *this;
}

template <const mpz_class& p>
Fp2Elem<p> Fp2Elem<p>::operator + (const Fp2Elem<p>& op) const
{
    return Fp2Elem(_real_ + op._real_, _imag_ + op._imag_);
}

template <const mpz_class& p>
Fp2Elem<p> Fp2Elem<p>::operator + (Fp2Elem<p>&& op) const
{
    return Fp2Elem(_real_ + op._real_, _imag_ + op._imag_);
}

template <const mpz_class& p>
Fp2Elem<p>& Fp2Elem<p>::operator += (const Fp2Elem<p>& op)
{
    _real_ = _real_ + op._real_;
    _imag_ = _imag_ + op._imag_;
    return *this;
}

template <const mpz_class& p>
Fp2Elem<p>& Fp2Elem<p>::operator += (Fp2Elem<p>&& op)
{
    _real_ = _real_ + op._real_;
    _imag_ = _imag_ + op._imag_;
    return *this;
}

template <const mpz_class& p>
Fp2Elem<p> Fp2Elem<p>::operator - (void) const
{
    return Fp2Elem(-_real_, -_imag_);
}

template <const mpz_class& p>
Fp2Elem<p> Fp2Elem<p>::operator - (const Fp2Elem<p>& op) const
{
    return Fp2Elem(_real_ - op._real_, _imag_ - op._imag_);
}

template <const mpz_class& p>
Fp2Elem<p> Fp2Elem<p>::operator - (Fp2Elem<p>&& op) const
{
    return Fp2Elem(_real_ - op._real_, _imag_ - op._imag_);
}

template <const mpz_class& p>
Fp2Elem<p>& Fp2Elem<p>::operator -= (const Fp2Elem<p>& op)
{
    _real_ = _real_ - op._real_;
    _imag_ = _imag_ - op._imag_;
    return *this;
}

template <const mpz_class& p>
Fp2Elem<p>& Fp2Elem<p>::operator -= (Fp2Elem<p>&& op)
{
    _real_ = _real_ - op._real_;
    _imag_ = _imag_ - op._imag_;
    return *this;
}

template <const mpz_class& p>
Fp2Elem<p> Fp2Elem<p>::operator * (const Fp2Elem<p>& op) const
{
    return Fp2Elem(_real_ * op._real_ - _imag_ * op._imag_, _real_ * op._imag_ + _imag_ * op._real_);
}
    
template <const mpz_class& p>
Fp2Elem<p> Fp2Elem<p>::operator * (Fp2Elem<p>&& op) const
{
    return Fp2Elem(_real_ * op._real_ - _imag_ * op._imag_, _real_ * op._imag_ + _imag_ * op._real_);
}

template <const mpz_class& p>
Fp2Elem<p>& Fp2Elem<p>::operator *= (const Fp2Elem<p>& op)
{
    Fp2Elem<p> OP(op);
    FpElem<p> real(_real_);
    _real_ = real * OP._real_ - _imag_ * OP._imag_;
    _imag_ = real * OP._imag_ + _imag_ * OP._real_;
    return *this;
}

template <const mpz_class& p>
Fp2Elem<p>& Fp2Elem<p>::operator *= (Fp2Elem<p>&& op)
{
    _real_ = _real_ * op._real_ - _imag_ * op._imag_;
    _imag_ = _real_ * op._imag_ + _imag_ * op._real_;
    return *this;
}

template <const mpz_class& p>
Fp2Elem<p> Fp2Elem<p>::invert(void) const
{
    FpElem<p> tmp = (_real_ * _real_ + _imag_ * _imag_).invert();
    Fp2Elem<p> res;
    res._real_ = _real_ * tmp;
    res._imag_ = - _imag_;
    res._imag_ *= tmp;
    return res; 
}

template <const mpz_class& p>
bool Fp2Elem<p>::operator == (Fp2Elem<p>& op)
{
    return (_real_ == op._real_) && (_imag_ == op._imag_);
}

template <const mpz_class& p>
std::string Fp2Elem<p>::get_real_str(void) const
{
    return _real_.get_str();
}

template <const mpz_class& p>
std::string Fp2Elem<p>::get_imag_str(void) const
{
    return _imag_.get_str();
}

template <const mpz_class& p>
std::string Fp2Elem<p>::get_p_str(void) const
{
    return p.get_str();
}

template <const mpz_class& p>
void Fp2Elem<p>::print(void) const
{
    std::cout << "(" <<_real_.get_str() << " + " << _imag_.get_str() << " * i) % " << p << std::endl;
}

template <const mpz_class& p>
std::string Fp2Elem<p>::get_str(void) const
{
    return "(" + _real_.get_str() + " + " + _imag_.get_str() + " * i) % " + p.get_str();
}
