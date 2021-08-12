/*-*- C++ -*-*/
#pragma once

namespace pre {

template <typename T, typename U, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator+=(
        ArrayView<T, R> lhs, ArrayView<U, R> rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs += *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator-=(
        ArrayView<T, R> lhs, ArrayView<U, R> rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs -= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator*=(
        ArrayView<T, R> lhs, ArrayView<U, R> rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs *= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator/=(
        ArrayView<T, R> lhs, ArrayView<U, R> rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs /= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator%=(
        ArrayView<T, R> lhs, ArrayView<U, R> rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs %= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator&=(
        ArrayView<T, R> lhs, ArrayView<U, R> rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs &= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator|=(
        ArrayView<T, R> lhs, ArrayView<U, R> rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs |= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator^=(
        ArrayView<T, R> lhs, ArrayView<U, R> rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs ^= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator>>=(
        ArrayView<T, R> lhs, ArrayView<U, R> rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs >>= *itrrhs;
    return lhs;
}

template <typename T, typename U, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator<<=(
        ArrayView<T, R> lhs, ArrayView<U, R> rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs <<= *itrrhs;
    return lhs;
}

template <typename T, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator+=(
        ArrayView<T, R> lhs, const std::type_identity_t<T>& rhs) noexcept {
    for (auto itr = lhs.begin(); itr != lhs.end(); ++itr)
        *itr += rhs;
    return lhs;
}

template <typename T, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator-=(
        ArrayView<T, R> lhs, const std::type_identity_t<T>& rhs) noexcept {
    for (auto itr = lhs.begin(); itr != lhs.end(); ++itr)
        *itr -= rhs;
    return lhs;
}

template <typename T, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator*=(
        ArrayView<T, R> lhs, const std::type_identity_t<T>& rhs) noexcept {
    for (auto itr = lhs.begin(); itr != lhs.end(); ++itr)
        *itr *= rhs;
    return lhs;
}

template <typename T, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator/=(
        ArrayView<T, R> lhs, const std::type_identity_t<T>& rhs) noexcept {
    for (auto itr = lhs.begin(); itr != lhs.end(); ++itr)
        *itr /= rhs;
    return lhs;
}

template <typename T, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator%=(
        ArrayView<T, R> lhs, const std::type_identity_t<T>& rhs) noexcept {
    for (auto itr = lhs.begin(); itr != lhs.end(); ++itr)
        *itr %= rhs;
    return lhs;
}

template <typename T, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator&=(
        ArrayView<T, R> lhs, const std::type_identity_t<T>& rhs) noexcept {
    for (auto itr = lhs.begin(); itr != lhs.end(); ++itr)
        *itr &= rhs;
    return lhs;
}

template <typename T, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator|=(
        ArrayView<T, R> lhs, const std::type_identity_t<T>& rhs) noexcept {
    for (auto itr = lhs.begin(); itr != lhs.end(); ++itr)
        *itr |= rhs;
    return lhs;
}

template <typename T, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator^=(
        ArrayView<T, R> lhs, const std::type_identity_t<T>& rhs) noexcept {
    for (auto itr = lhs.begin(); itr != lhs.end(); ++itr)
        *itr ^= rhs;
    return lhs;
}

template <typename T, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator>>=(
        ArrayView<T, R> lhs, const std::type_identity_t<T>& rhs) noexcept {
    for (auto itr = lhs.begin(); itr != lhs.end(); ++itr)
        *itr >>= rhs;
    return lhs;
}

template <typename T, size_t R>
[[gnu::always_inline]] constexpr ArrayView<T, R> operator<<=(
        ArrayView<T, R> lhs, const std::type_identity_t<T>& rhs) noexcept {
    for (auto itr = lhs.begin(); itr != lhs.end(); ++itr)
        *itr <<= rhs;
    return lhs;
}

} // namespace pre
