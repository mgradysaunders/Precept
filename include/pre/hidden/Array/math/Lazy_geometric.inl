/*-*- C++ -*-*/
#pragma once

namespace pre {

/// Euclidean length.
template <typename Func, size_t Rank>
[[gnu::always_inline]] inline auto length(
    const LazyArray<Func, Rank>& arr) noexcept {
    return pre::sqrt(pre::sum(pre::norm(arr)));
}

/// Euclidean length squared.
template <typename Func, size_t Rank>
[[gnu::always_inline]] inline auto length2(
    const LazyArray<Func, Rank>& arr) noexcept {
    return pre::sum(pre::norm(arr));
}

/// Normalize by Euclidean length.
template <typename Func>
[[gnu::always_inline]] inline auto normalize(
    const LazyArray<Func, 1>& arr) noexcept {
    using Result = typename LazyArray<Func, 1>::Result;
    using Float = to_floating_point_t<Result>;
    Float fac = 1 / pre::length(arr);
    if (not std::isfinite(fac))
        fac = 0;
    return LazyArray(
        [&, fac](auto ind) noexcept { return fac * arr(ind); }, arr.sizes);
}

} // namespace pre
