/*-*- C++ -*-*/
#pragma once

namespace pre {

template <typename F, size_t R>
[[gnu::always_inline]] inline auto abs(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::abs(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto arg(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::arg(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto real(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::real(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto imag(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::imag(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto conj(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::conj(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto norm(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::norm(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto fabs(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::fabs(arrx(k)); }, arrx.sizes);
}

template <typename F, typename G, typename H, size_t R>
[[gnu::always_inline]] inline auto fma(
        const LazyArray<F, R>& arrx,
        const LazyArray<G, R>& arry,
        const LazyArray<H, R>& arrz) noexcept {
    return LazyArray(
            [&](auto k) { return pre::fma(arrx(k), arry(k), arrz(k)); },
            arrx.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] inline auto fmin(
        const LazyArray<F, R>& arrx, const LazyArray<G, R>& arry) noexcept {
    return LazyArray(
            [&](auto k) { return pre::fmin(arrx(k), arry(k)); }, arrx.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] inline auto fmax(
        const LazyArray<F, R>& arrx, const LazyArray<G, R>& arry) noexcept {
    return LazyArray(
            [&](auto k) { return pre::fmax(arrx(k), arry(k)); }, arrx.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] inline auto fdim(
        const LazyArray<F, R>& arrx, const LazyArray<G, R>& arry) noexcept {
    return LazyArray(
            [&](auto k) { return pre::fdim(arrx(k), arry(k)); }, arrx.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] inline auto fmod(
        const LazyArray<F, R>& arrx, const LazyArray<G, R>& arry) noexcept {
    return LazyArray(
            [&](auto k) { return pre::fmod(arrx(k), arry(k)); }, arrx.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] inline auto remainder(
        const LazyArray<F, R>& arrx, const LazyArray<G, R>& arry) noexcept {
    return LazyArray(
            [&](auto k) { return pre::remainder(arrx(k), arry(k)); },
            arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto nearbyint(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray(
            [&](auto k) { return pre::nearbyint(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto floor(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::floor(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto ceil(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::ceil(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto trunc(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::trunc(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto round(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::round(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto rint(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::rint(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto lrint(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::lrint(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto llrint(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::llrint(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto lround(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::lround(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto llround(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray(
            [&](auto k) { return pre::llround(arrx(k)); }, arrx.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] inline auto nextafter(
        const LazyArray<F, R>& arrx, const LazyArray<G, R>& arry) noexcept {
    return LazyArray(
            [&](auto k) { return pre::nextafter(arrx(k), arry(k)); },
            arrx.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] inline auto nexttoward(
        const LazyArray<F, R>& arrx, const LazyArray<G, R>& arry) noexcept {
    return LazyArray(
            [&](auto k) { return pre::nexttoward(arrx(k), arry(k)); },
            arrx.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] inline auto copysign(
        const LazyArray<F, R>& arrx, const LazyArray<G, R>& arry) noexcept {
    return LazyArray(
            [&](auto k) { return pre::copysign(arrx(k), arry(k)); },
            arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto exp(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::exp(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto log(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::log(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto exp2(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::exp2(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto log2(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::log2(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto log10(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::log10(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto expm1(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::expm1(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto log1p(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::log1p(arrx(k)); }, arrx.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] inline auto pow(
        const LazyArray<F, R>& arrx, const LazyArray<G, R>& arry) noexcept {
    return LazyArray(
            [&](auto k) { return pre::pow(arrx(k), arry(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto sqrt(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::sqrt(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto cbrt(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::cbrt(arrx(k)); }, arrx.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] inline auto hypot(
        const LazyArray<F, R>& arrx, const LazyArray<G, R>& arry) noexcept {
    return LazyArray(
            [&](auto k) { return pre::hypot(arrx(k), arry(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto erf(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::erf(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto erfc(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::erfc(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto lgamma(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::lgamma(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto tgamma(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::tgamma(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto sin(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::sin(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto cos(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::cos(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto tan(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::tan(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto asin(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::asin(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto acos(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::acos(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto atan(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::atan(arrx(k)); }, arrx.sizes);
}

template <typename G, typename F, size_t R>
[[gnu::always_inline]] inline auto atan2(
        const LazyArray<G, R>& arry, const LazyArray<F, R>& arrx) noexcept {
    return LazyArray(
            [&](auto k) { return pre::atan2(arry(k), arrx(k)); }, arry.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto sinh(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::sinh(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto cosh(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::cosh(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto tanh(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::tanh(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto asinh(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::asinh(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto acosh(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::acosh(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto atanh(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::atanh(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto csc(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::csc(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto sec(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::sec(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto cot(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::cot(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto csch(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::csch(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto sech(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::sech(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto coth(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::coth(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto acsc(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::acsc(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto asec(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::asec(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto acot(const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::acot(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto acsch(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::acsch(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto asech(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::asech(arrx(k)); }, arrx.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] inline auto acoth(
        const LazyArray<F, R>& arrx) noexcept {
    return LazyArray([&](auto k) { return pre::acoth(arrx(k)); }, arrx.sizes);
}

} // namespace pre
