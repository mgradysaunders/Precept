/*-*- C++ -*-*/
#pragma once

namespace pre {

/// Dual number implementation of `pre::exp()`.
template <typename T>
inline Dual<T> exp(const Dual<T>& x) {
    return {pre::exp(x.real()), pre::exp(x.real()) * x.dual()};
}

/// Dual number implementation of `pre::log()`.
template <typename T>
inline Dual<T> log(const Dual<T>& x) {
    return {pre::log(x.real()), x.dual() / x.real()};
}

/// Dual number implementation of `pre::exp2()`.
template <typename T>
inline Dual<T> exp2(const Dual<T>& x) {
    return {pre::exp2(x.real()), pre::numeric_constants<T>::M_ln2() *
                                         pre::exp2(x.real()) * x.dual()};
}

/// Dual number implementation of `pre::log2()`.
template <typename T>
inline Dual<T> log2(const Dual<T>& x) {
    return {pre::log2(x.real()),
            x.dual() / (pre::numeric_constants<T>::M_ln2() * x.real())};
}

/// Dual number implementation of `pre::log10()`.
template <typename T>
inline Dual<T> log10(const Dual<T>& x) {
    return {pre::log10(x.real()),
            x.dual() / (pre::numeric_constants<T>::M_ln10() * x.real())};
}

/// Dual number implementation of `pre::expm1()`.
template <typename T>
inline Dual<T> expm1(const Dual<T>& x) {
    return {pre::expm1(x.real()), pre::exp(x.real()) * x.dual()};
}

/// Dual number implementation of `pre::log1p()`.
template <typename T>
inline Dual<T> log1p(const Dual<T>& x) {
    return {pre::log1p(x.real()), x.dual() / (T(1) + x.real())};
}

/// Dual number implementation of `pre::sqrt()`.
template <typename T>
inline Dual<T> sqrt(const Dual<T>& x) {
    return {pre::sqrt(x.real()), x.dual() / (T(2) * pre::sqrt(x.real()))};
}

/// Dual number implementation of `pre::cbrt()`.
template <typename T>
inline Dual<T> cbrt(const Dual<T>& x) {
    return {pre::cbrt(x.real()),
            x.dual() / (T(3) * pre::nthpow(pre::cbrt(x.real()), 2))};
}

/// Dual number implementation of `pre::erf()`.
template <typename T>
inline Dual<T> erf(const Dual<T>& x) {
    return {pre::erf(x.real()), pre::numeric_constants<T>::M_2_sqrtpi() *
                                        pre::exp(-pre::nthpow(x.real(), 2)) *
                                        x.dual()};
}

/// Dual number implementation of `pre::erfc()`.
template <typename T>
inline Dual<T> erfc(const Dual<T>& x) {
    return {pre::erfc(x.real()), -pre::numeric_constants<T>::M_2_sqrtpi() *
                                         pre::exp(-pre::nthpow(x.real(), 2)) *
                                         x.dual()};
}

/// Dual number implementation of `pre::sin()`.
template <typename T>
inline Dual<T> sin(const Dual<T>& x) {
    return {pre::sin(x.real()), pre::cos(x.real()) * x.dual()};
}

/// Dual number implementation of `pre::cos()`.
template <typename T>
inline Dual<T> cos(const Dual<T>& x) {
    return {pre::cos(x.real()), -pre::sin(x.real()) * x.dual()};
}

/// Dual number implementation of `pre::tan()`.
template <typename T>
inline Dual<T> tan(const Dual<T>& x) {
    return {pre::tan(x.real()), x.dual() / pre::nthpow(pre::cos(x.real()), 2)};
}

/// Dual number implementation of `pre::asin()`.
template <typename T>
inline Dual<T> asin(const Dual<T>& x) {
    return {pre::asin(x.real()),
            x.dual() / pre::sqrt(T(1) - pre::nthpow(x.real(), 2))};
}

/// Dual number implementation of `pre::acos()`.
template <typename T>
inline Dual<T> acos(const Dual<T>& x) {
    return {pre::acos(x.real()),
            x.dual() / -pre::sqrt(T(1) - pre::nthpow(x.real(), 2))};
}

/// Dual number implementation of `pre::atan()`.
template <typename T>
inline Dual<T> atan(const Dual<T>& x) {
    return {pre::atan(x.real()), x.dual() / (T(1) + pre::nthpow(x.real(), 2))};
}

/// Dual number implementation of `pre::sinh()`.
template <typename T>
inline Dual<T> sinh(const Dual<T>& x) {
    return {pre::sinh(x.real()), pre::cosh(x.real()) * x.dual()};
}

/// Dual number implementation of `pre::cosh()`.
template <typename T>
inline Dual<T> cosh(const Dual<T>& x) {
    return {pre::cosh(x.real()), pre::sinh(x.real()) * x.dual()};
}

/// Dual number implementation of `pre::tanh()`.
template <typename T>
inline Dual<T> tanh(const Dual<T>& x) {
    return {pre::tanh(x.real()),
            x.dual() / pre::nthpow(pre::cosh(x.real()), 2)};
}

/// Dual number implementation of `pre::asinh()`.
template <typename T>
inline Dual<T> asinh(const Dual<T>& x) {
    return {pre::asinh(x.real()),
            x.dual() / pre::sqrt(pre::nthpow(x.real(), 2) + T(1))};
}

/// Dual number implementation of `pre::acosh()`.
template <typename T>
inline Dual<T> acosh(const Dual<T>& x) {
    return {pre::acosh(x.real()),
            x.dual() / pre::sqrt(pre::nthpow(x.real(), 2) - T(1))};
}

/// Dual number implementation of `pre::atanh()`.
template <typename T>
inline Dual<T> atanh(const Dual<T>& x) {
    return {pre::atanh(x.real()),
            x.dual() / (T(1) - pre::nthpow(x.real(), 2))};
}

/// Dual number implementation of `pre::csc()`.
template <typename T>
[[gnu::always_inline]] inline Dual<T> csc(const Dual<T>& x) {
    return T(1) / pre::sin(x);
}

/// Dual number implementation of `pre::sec()`.
template <typename T>
[[gnu::always_inline]] inline Dual<T> sec(const Dual<T>& x) {
    return T(1) / pre::cos(x);
}

/// Dual number implementation of `pre::cot()`.
template <typename T>
[[gnu::always_inline]] inline Dual<T> cot(const Dual<T>& x) {
    return T(1) / pre::tan(x);
}

/// Dual number implementation of `pre::csch()`.
template <typename T>
[[gnu::always_inline]] inline Dual<T> csch(const Dual<T>& x) {
    return T(1) / pre::sinh(x);
}

/// Dual number implementation of `pre::sech()`.
template <typename T>
[[gnu::always_inline]] inline Dual<T> sech(const Dual<T>& x) {
    return T(1) / pre::cosh(x);
}

/// Dual number implementation of `pre::coth()`.
template <typename T>
[[gnu::always_inline]] inline Dual<T> coth(const Dual<T>& x) {
    return T(1) / pre::tanh(x);
}

/// Dual number implementation of `pre::acsc()`.
template <typename T>
[[gnu::always_inline]] inline Dual<T> acsc(const Dual<T>& x) {
    return pre::asin(T(1) / x);
}

/// Dual number implementation of `pre::asec()`.
template <typename T>
[[gnu::always_inline]] inline Dual<T> asec(const Dual<T>& x) {
    return pre::acos(T(1) / x);
}

/// Dual number implementation of `pre::acot()`.
template <typename T>
[[gnu::always_inline]] inline Dual<T> acot(const Dual<T>& x) {
    return pre::atan(T(1) / x);
}

/// Dual number implementation of `pre::acsch()`.
template <typename T>
[[gnu::always_inline]] inline Dual<T> acsch(const Dual<T>& x) {
    return pre::asinh(T(1) / x);
}

/// Dual number implementation of `pre::asech()`.
template <typename T>
[[gnu::always_inline]] inline Dual<T> asech(const Dual<T>& x) {
    return pre::acosh(T(1) / x);
}

/// Dual number implementation of `pre::acoth()`.
template <typename T>
[[gnu::always_inline]] inline Dual<T> acoth(const Dual<T>& x) {
    return pre::atanh(T(1) / x);
}

} // namespace pre
