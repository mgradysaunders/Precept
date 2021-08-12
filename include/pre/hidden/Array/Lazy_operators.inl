/*-*- C++ -*-*/
#pragma once

namespace pre {

template <typename F, size_t R>
[[gnu::always_inline]] constexpr auto operator+(
        const LazyArray<F, R>& arr) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return +arr(k); }, arr.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] constexpr auto operator-(
        const LazyArray<F, R>& arr) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return -arr(k); }, arr.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] constexpr auto operator~(
        const LazyArray<F, R>& arr) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return ~arr(k); }, arr.sizes);
}

template <typename F, size_t R>
[[gnu::always_inline]] constexpr auto operator!(
        const LazyArray<F, R>& arr) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return !arr(k); }, arr.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] constexpr auto operator+(
        const LazyArray<F, R>& lhs, const LazyArray<G, R>& rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) + rhs(k); },
            lhs.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] constexpr auto operator-(
        const LazyArray<F, R>& lhs, const LazyArray<G, R>& rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) - rhs(k); },
            lhs.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] constexpr auto operator*(
        const LazyArray<F, R>& lhs, const LazyArray<G, R>& rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) * rhs(k); },
            lhs.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] constexpr auto operator/(
        const LazyArray<F, R>& lhs, const LazyArray<G, R>& rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) / rhs(k); },
            lhs.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] constexpr auto operator%(
        const LazyArray<F, R>& lhs, const LazyArray<G, R>& rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) % rhs(k); },
            lhs.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] constexpr auto operator&(
        const LazyArray<F, R>& lhs, const LazyArray<G, R>& rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) & rhs(k); },
            lhs.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] constexpr auto operator|(
        const LazyArray<F, R>& lhs, const LazyArray<G, R>& rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) | rhs(k); },
            lhs.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] constexpr auto operator^(
        const LazyArray<F, R>& lhs, const LazyArray<G, R>& rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) ^ rhs(k); },
            lhs.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] constexpr auto operator>>(
        const LazyArray<F, R>& lhs, const LazyArray<G, R>& rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) >> rhs(k); },
            lhs.sizes);
}

template <typename F, typename G, size_t R>
[[gnu::always_inline]] constexpr auto operator<<(
        const LazyArray<F, R>& lhs, const LazyArray<G, R>& rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) << rhs(k); },
            lhs.sizes);
}

template <typename F, concepts::arithmetic_or_complex T, size_t R>
[[gnu::always_inline]] constexpr auto operator+(
        const LazyArray<F, R>& lhs, const T& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) + rhs; },
            lhs.sizes);
}

template <typename F, concepts::arithmetic_or_complex T, size_t R>
[[gnu::always_inline]] constexpr auto operator-(
        const LazyArray<F, R>& lhs, const T& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) - rhs; },
            lhs.sizes);
}

template <typename F, concepts::arithmetic_or_complex T, size_t R>
[[gnu::always_inline]] constexpr auto operator*(
        const LazyArray<F, R>& lhs, const T& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) * rhs; },
            lhs.sizes);
}

template <typename F, concepts::arithmetic_or_complex T, size_t R>
[[gnu::always_inline]] constexpr auto operator/(
        const LazyArray<F, R>& lhs, const T& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) / rhs; },
            lhs.sizes);
}

template <typename F, concepts::arithmetic_or_complex T, size_t R>
[[gnu::always_inline]] constexpr auto operator%(
        const LazyArray<F, R>& lhs, const T& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) % rhs; },
            lhs.sizes);
}

template <typename F, concepts::arithmetic_or_complex T, size_t R>
[[gnu::always_inline]] constexpr auto operator&(
        const LazyArray<F, R>& lhs, const T& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) & rhs; },
            lhs.sizes);
}

template <typename F, concepts::arithmetic_or_complex T, size_t R>
[[gnu::always_inline]] constexpr auto operator|(
        const LazyArray<F, R>& lhs, const T& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) | rhs; },
            lhs.sizes);
}

template <typename F, concepts::arithmetic_or_complex T, size_t R>
[[gnu::always_inline]] constexpr auto operator^(
        const LazyArray<F, R>& lhs, const T& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) ^ rhs; },
            lhs.sizes);
}

template <typename F, concepts::arithmetic_or_complex T, size_t R>
[[gnu::always_inline]] constexpr auto operator>>(
        const LazyArray<F, R>& lhs, const T& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) >> rhs; },
            lhs.sizes);
}

template <typename F, concepts::arithmetic_or_complex T, size_t R>
[[gnu::always_inline]] constexpr auto operator<<(
        const LazyArray<F, R>& lhs, const T& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs(k) << rhs; },
            lhs.sizes);
}

template <concepts::arithmetic_or_complex T, typename F, size_t R>
[[gnu::always_inline]] constexpr auto operator+(
        const T& lhs, const LazyArray<F, R>& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs + rhs(k); },
            rhs.sizes);
}

template <concepts::arithmetic_or_complex T, typename F, size_t R>
[[gnu::always_inline]] constexpr auto operator-(
        const T& lhs, const LazyArray<F, R>& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs - rhs(k); },
            rhs.sizes);
}

template <concepts::arithmetic_or_complex T, typename F, size_t R>
[[gnu::always_inline]] constexpr auto operator*(
        const T& lhs, const LazyArray<F, R>& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs * rhs(k); },
            rhs.sizes);
}

template <concepts::arithmetic_or_complex T, typename F, size_t R>
[[gnu::always_inline]] constexpr auto operator/(
        const T& lhs, const LazyArray<F, R>& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs / rhs(k); },
            rhs.sizes);
}

template <concepts::arithmetic_or_complex T, typename F, size_t R>
[[gnu::always_inline]] constexpr auto operator%(
        const T& lhs, const LazyArray<F, R>& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs % rhs(k); },
            rhs.sizes);
}

template <concepts::arithmetic_or_complex T, typename F, size_t R>
[[gnu::always_inline]] constexpr auto operator&(
        const T& lhs, const LazyArray<F, R>& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs & rhs(k); },
            rhs.sizes);
}

template <concepts::arithmetic_or_complex T, typename F, size_t R>
[[gnu::always_inline]] constexpr auto operator|(
        const T& lhs, const LazyArray<F, R>& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs | rhs(k); },
            rhs.sizes);
}

template <concepts::arithmetic_or_complex T, typename F, size_t R>
[[gnu::always_inline]] constexpr auto operator^(
        const T& lhs, const LazyArray<F, R>& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs ^ rhs(k); },
            rhs.sizes);
}

template <concepts::arithmetic_or_complex T, typename F, size_t R>
[[gnu::always_inline]] constexpr auto operator>>(
        const T& lhs, const LazyArray<F, R>& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs >> rhs(k); },
            rhs.sizes);
}

template <concepts::arithmetic_or_complex T, typename F, size_t R>
[[gnu::always_inline]] constexpr auto operator<<(
        const T& lhs, const LazyArray<F, R>& rhs) noexcept {
    return LazyArray(
            [&](auto k) constexpr noexcept { return lhs << rhs(k); },
            rhs.sizes);
}

} // namespace pre
