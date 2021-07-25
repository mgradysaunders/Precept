/*-*- C++ -*-*/
#pragma once

namespace pre {

template <typename T, size_t... N>
[[gnu::always_inline]] constexpr auto operator+(
        const Array<T, N...>& arr) noexcept {
    using U = decltype(+T());
    Array<U, N...> res;
    auto itrarr = arr.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrarr, ++itrres)
        *itrres = +*itrarr;
    return res;
}

template <typename T, size_t... N>
[[gnu::always_inline]] constexpr auto operator-(
        const Array<T, N...>& arr) noexcept {
    using U = decltype(-T());
    Array<U, N...> res;
    auto itrarr = arr.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrarr, ++itrres)
        *itrres = -*itrarr;
    return res;
}

template <typename T, size_t... N>
[[gnu::always_inline]] constexpr auto operator~(
        const Array<T, N...>& arr) noexcept {
    using U = decltype(~T());
    Array<U, N...> res;
    auto itrarr = arr.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrarr, ++itrres)
        *itrres = ~*itrarr;
    return res;
}

template <typename T, size_t... N>
[[gnu::always_inline]] constexpr auto operator!(
        const Array<T, N...>& arr) noexcept {
    using U = decltype(!T());
    Array<U, N...> res;
    auto itrarr = arr.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrarr, ++itrres)
        *itrres = !*itrarr;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator+(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() + U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs + *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator-(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() - U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs - *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator*(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() * U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs * *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator/(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() / U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs / *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator%(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() % U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs % *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator&(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() & U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs & *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator|(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() | U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs | *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator^(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() ^ U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs ^ *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator>>(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() >> U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs >> *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator<<(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() << U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs << *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator==(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() == U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs == *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator!=(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() != U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs != *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator<(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() < U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs < *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator>(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() > U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs > *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator<=(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() <= U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs <= *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator>=(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() >= U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs >= *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator&&(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() && U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs && *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator||(
        const Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() || U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs || *itrrhs;
    return res;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator+=(
        Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs += *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator-=(
        Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs -= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator*=(
        Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs *= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator/=(
        Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs /= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator%=(
        Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs %= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator&=(
        Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs &= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator|=(
        Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs |= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator^=(
        Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs ^= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator>>=(
        Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs >>= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator<<=(
        Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept {
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs <<= *itrrhs;
    return lhs;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator+(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() + U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs + rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator-(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() - U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs - rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator*(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() * U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs * rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator/(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() / U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs / rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator%(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() % U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs % rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator&(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() & U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs & rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator|(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() | U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs | rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator^(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() ^ U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs ^ rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator>>(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() >> U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs >> rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator<<(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() << U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs << rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator==(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() == U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs == rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator!=(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() != U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs != rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator<(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() < U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs < rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator>(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() > U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs > rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator<=(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() <= U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs <= rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator>=(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() >= U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs >= rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator&&(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() && U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs && rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr auto operator||(
        const Array<T, N...>& lhs, const U& rhs) noexcept {
    using V = decltype(T() || U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs || rhs;
    return res;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator+=(
        Array<T, N...>& lhs, const U& rhs) noexcept {
    for (auto itrlhs = lhs.begin(); itrlhs != lhs.end(); ++itrlhs)
        *itrlhs += rhs;
    return lhs;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator-=(
        Array<T, N...>& lhs, const U& rhs) noexcept {
    for (auto itrlhs = lhs.begin(); itrlhs != lhs.end(); ++itrlhs)
        *itrlhs -= rhs;
    return lhs;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator*=(
        Array<T, N...>& lhs, const U& rhs) noexcept {
    for (auto itrlhs = lhs.begin(); itrlhs != lhs.end(); ++itrlhs)
        *itrlhs *= rhs;
    return lhs;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator/=(
        Array<T, N...>& lhs, const U& rhs) noexcept {
    for (auto itrlhs = lhs.begin(); itrlhs != lhs.end(); ++itrlhs)
        *itrlhs /= rhs;
    return lhs;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator%=(
        Array<T, N...>& lhs, const U& rhs) noexcept {
    for (auto itrlhs = lhs.begin(); itrlhs != lhs.end(); ++itrlhs)
        *itrlhs %= rhs;
    return lhs;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator&=(
        Array<T, N...>& lhs, const U& rhs) noexcept {
    for (auto itrlhs = lhs.begin(); itrlhs != lhs.end(); ++itrlhs)
        *itrlhs &= rhs;
    return lhs;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator|=(
        Array<T, N...>& lhs, const U& rhs) noexcept {
    for (auto itrlhs = lhs.begin(); itrlhs != lhs.end(); ++itrlhs)
        *itrlhs |= rhs;
    return lhs;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator^=(
        Array<T, N...>& lhs, const U& rhs) noexcept {
    for (auto itrlhs = lhs.begin(); itrlhs != lhs.end(); ++itrlhs)
        *itrlhs ^= rhs;
    return lhs;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator>>=(
        Array<T, N...>& lhs, const U& rhs) noexcept {
    for (auto itrlhs = lhs.begin(); itrlhs != lhs.end(); ++itrlhs)
        *itrlhs >>= rhs;
    return lhs;
}

template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]] constexpr Array<T, N...>& operator<<=(
        Array<T, N...>& lhs, const U& rhs) noexcept {
    for (auto itrlhs = lhs.begin(); itrlhs != lhs.end(); ++itrlhs)
        *itrlhs <<= rhs;
    return lhs;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator+(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() + U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs + *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator-(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() - U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs - *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator*(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() * U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs * *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator/(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() / U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs / *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator%(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() % U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs % *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator&(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() & U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs & *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator|(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() | U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs | *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator^(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() ^ U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs ^ *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator>>(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() >> U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs >> *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator<<(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() << U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs << *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator==(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() == U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs == *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator!=(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() != U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs != *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator<(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() < U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs < *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator>(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() > U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs > *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator<=(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() <= U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs <= *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator>=(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() >= U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs >= *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator&&(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() && U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs && *itrrhs;
    return res;
}

template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]] constexpr auto operator||(
        const T& lhs, const Array<U, N...>& rhs) noexcept {
    using V = decltype(T() || U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs || *itrrhs;
    return res;
}

} // namespace pre
