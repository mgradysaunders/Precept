/*-*- C++ -*-*/
#pragma once

namespace pre {

template <typename T, size_t... N>
inline auto abs(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::abs(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::abs(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto arg(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::arg(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::arg(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto real(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::real(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::real(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto imag(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::imag(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::imag(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto conj(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::conj(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::conj(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto norm(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::norm(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::norm(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto fabs(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::fabs(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::fabs(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto fma(
        const Array<T, N...>& arrx,
        const Array<T, N...>& arry,
        const Array<T, N...>& arrz) noexcept {
    Array<decltype(pre::fma(T(), T(), T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarry = arry.begin();
    auto itrarrz = arrz.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarry, ++itrarrz, ++itrres)
        *itrres = pre::fma(*itrarrx, *itrarry, *itrarrz);
    return res;
}

template <typename T, size_t... N>
inline auto fmin(
        const Array<T, N...>& arrx, const Array<T, N...>& arry) noexcept {
    Array<decltype(pre::fmin(T(), T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarry = arry.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarry, ++itrres)
        *itrres = pre::fmin(*itrarrx, *itrarry);
    return res;
}

template <typename T, size_t... N>
inline auto fmax(
        const Array<T, N...>& arrx, const Array<T, N...>& arry) noexcept {
    Array<decltype(pre::fmax(T(), T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarry = arry.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarry, ++itrres)
        *itrres = pre::fmax(*itrarrx, *itrarry);
    return res;
}

template <typename T, size_t... N>
inline auto fdim(
        const Array<T, N...>& arrx, const Array<T, N...>& arry) noexcept {
    Array<decltype(pre::fdim(T(), T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarry = arry.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarry, ++itrres)
        *itrres = pre::fdim(*itrarrx, *itrarry);
    return res;
}

template <typename T, size_t... N>
inline auto fmod(
        const Array<T, N...>& arrx, const Array<T, N...>& arry) noexcept {
    Array<decltype(pre::fmod(T(), T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarry = arry.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarry, ++itrres)
        *itrres = pre::fmod(*itrarrx, *itrarry);
    return res;
}

template <typename T, size_t... N>
inline auto remquo(
        const Array<T, N...>& arrx,
        const Array<T, N...>& arry,
        Array<int, N...>* arrq) noexcept {
    Array<decltype(pre::remquo(T(), T(), std::declval<int*>())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarry = arry.begin();
    auto itrarrq = arrq->begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarry, ++itrarrq, ++itrres)
        *itrres = pre::remquo(*itrarrx, *itrarry, &(*itrarrq));
    return res;
}

template <typename T, size_t... N>
inline auto remainder(
        const Array<T, N...>& arrx, const Array<T, N...>& arry) noexcept {
    Array<decltype(pre::remainder(T(), T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarry = arry.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarry, ++itrres)
        *itrres = pre::remainder(*itrarrx, *itrarry);
    return res;
}

template <typename T, size_t... N>
inline auto nearbyint(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::nearbyint(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::nearbyint(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto floor(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::floor(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::floor(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto ceil(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::ceil(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::ceil(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto trunc(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::trunc(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::trunc(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto round(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::round(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::round(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto rint(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::rint(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::rint(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto lrint(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::lrint(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::lrint(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto llrint(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::llrint(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::llrint(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto lround(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::lround(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::lround(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto llround(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::llround(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::llround(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto frexp(
        const Array<T, N...>& arrx, Array<int, N...>* arrp) noexcept {
    Array<decltype(pre::frexp(T(), std::declval<int*>())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarrp = arrp->begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarrp, ++itrres)
        *itrres = pre::frexp(*itrarrx, &(*itrarrp));
    return res;
}

template <typename T, size_t... N>
inline auto ldexp(
        const Array<T, N...>& arrx, const Array<int, N...>& arrp) noexcept {
    Array<decltype(pre::ldexp(T(), std::declval<int>())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarrp = arrp.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarrp, ++itrres)
        *itrres = pre::ldexp(*itrarrx, *itrarrp);
    return res;
}

template <typename T, size_t... N>
inline auto logb(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::logb(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::logb(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto ilogb(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::ilogb(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::ilogb(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto scalbn(
        const Array<T, N...>& arrx, const Array<int, N...>& arrp) noexcept {
    Array<decltype(pre::scalbn(T(), std::declval<int>())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarrp = arrp.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarrp, ++itrres)
        *itrres = pre::scalbn(*itrarrx, *itrarrp);
    return res;
}

template <typename T, size_t... N>
inline auto scalbln(
        const Array<T, N...>& arrx, const Array<long, N...>& arrp) noexcept {
    Array<decltype(pre::scalbln(T(), std::declval<long>())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarrp = arrp.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarrp, ++itrres)
        *itrres = pre::scalbln(*itrarrx, *itrarrp);
    return res;
}

template <typename T, size_t... N>
inline auto modf(const Array<T, N...>& arrx, Array<T, N...>* arrp) noexcept {
    Array<decltype(pre::modf(T(), std::declval<T*>())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarrp = arrp->begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarrp, ++itrres)
        *itrres = pre::modf(*itrarrx, &(*itrarrp));
    return res;
}

template <typename T, size_t... N>
inline auto nextafter(
        const Array<T, N...>& arrx, const Array<T, N...>& arry) noexcept {
    Array<decltype(pre::nextafter(T(), T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarry = arry.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarry, ++itrres)
        *itrres = pre::nextafter(*itrarrx, *itrarry);
    return res;
}

template <typename T, size_t... N>
inline auto nexttoward(
        const Array<T, N...>& arrx,
        const Array<long double, N...>& arry) noexcept {
    Array<decltype(pre::nexttoward(T(), std::declval<long double>())), N...>
            res;
    auto itrarrx = arrx.begin();
    auto itrarry = arry.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarry, ++itrres)
        *itrres = pre::nexttoward(*itrarrx, *itrarry);
    return res;
}

template <typename T, size_t... N>
inline auto copysign(
        const Array<T, N...>& arrx, const Array<T, N...>& arry) noexcept {
    Array<decltype(pre::copysign(T(), T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarry = arry.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarry, ++itrres)
        *itrres = pre::copysign(*itrarrx, *itrarry);
    return res;
}

template <typename T, size_t... N>
inline auto signbit(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::signbit(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::signbit(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto isnan(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::isnan(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::isnan(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto isinf(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::isinf(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::isinf(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto isfinite(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::isfinite(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::isfinite(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto isnormal(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::isnormal(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::isnormal(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto exp(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::exp(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::exp(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto log(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::log(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::log(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto exp2(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::exp2(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::exp2(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto log2(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::log2(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::log2(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto log10(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::log10(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::log10(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto expm1(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::expm1(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::expm1(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto log1p(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::log1p(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::log1p(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto pow(
        const Array<T, N...>& arrx, const Array<T, N...>& arry) noexcept {
    Array<decltype(pre::pow(T(), T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarry = arry.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarry, ++itrres)
        *itrres = pre::pow(*itrarrx, *itrarry);
    return res;
}

template <typename T, size_t... N>
inline auto sqrt(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::sqrt(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::sqrt(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto cbrt(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::cbrt(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::cbrt(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto hypot(
        const Array<T, N...>& arrx, const Array<T, N...>& arry) noexcept {
    Array<decltype(pre::hypot(T(), T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarry = arry.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarry, ++itrres)
        *itrres = pre::hypot(*itrarrx, *itrarry);
    return res;
}

template <typename T, size_t... N>
inline auto erf(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::erf(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::erf(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto erfc(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::erfc(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::erfc(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto lgamma(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::lgamma(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::lgamma(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto tgamma(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::tgamma(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::tgamma(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto sin(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::sin(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::sin(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto cos(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::cos(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::cos(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto tan(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::tan(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::tan(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto asin(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::asin(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::asin(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto acos(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::acos(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::acos(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto atan(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::atan(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::atan(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto atan2(
        const Array<T, N...>& arry, const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::atan2(T(), T())), N...> res;
    auto itrarry = arry.begin();
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarry, ++itrarrx, ++itrres)
        *itrres = pre::atan2(*itrarry, *itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto sinh(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::sinh(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::sinh(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto cosh(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::cosh(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::cosh(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto tanh(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::tanh(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::tanh(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto asinh(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::asinh(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::asinh(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto acosh(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::acosh(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::acosh(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto atanh(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::atanh(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::atanh(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto csc(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::csc(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::csc(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto sec(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::sec(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::sec(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto cot(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::cot(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::cot(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto csch(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::csch(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::csch(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto sech(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::sech(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::sech(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto coth(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::coth(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::coth(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto acsc(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::acsc(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::acsc(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto asec(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::asec(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::asec(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto acot(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::acot(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::acot(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto acsch(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::acsch(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::acsch(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto asech(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::asech(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::asech(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto acoth(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::acoth(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::acoth(*itrarrx);
    return res;
}

} // namespace pre
