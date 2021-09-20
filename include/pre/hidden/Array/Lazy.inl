/*-*- C++ -*-*/
#pragma once

namespace pre {

template <typename Func, size_t Rank>
struct LazyArray {
    using Result = std::decay_t<std::invoke_result_t<Func, ArrayIndex<Rank>>>;
    Func func;
    ArrayIndex<Rank> sizes;
    explicit constexpr LazyArray(Func&& f, const ArrayIndex<Rank>& s) noexcept
        : func(f), sizes(s) {
    }
    template <typename... Args>
    constexpr auto operator()(Args&&... args) const noexcept {
        return std::invoke(
            func, ArrayIndex<Rank>{std::forward<Args>(args)...});
    }
};

template <typename Func, size_t Rank>
LazyArray(Func&&, const ArrayIndex<Rank>&) -> LazyArray<Func, Rank>;

/// Find minimum.
template <typename Func>
[[gnu::always_inline]] constexpr auto min(
    const LazyArray<Func, 1>& arr) noexcept {
    using Result = typename LazyArray<Func, 1>::Result;
    if (arr.sizes[0] == 0)
        return Result{};
    Result res = arr(0);
    for (ssize_t ind = 1; ind < arr.sizes[0]; ind++)
        res = std::min(res, arr(ind));
    return res;
}

/// Find maximum.
template <typename Func>
[[gnu::always_inline]] constexpr auto max(
    const LazyArray<Func, 1>& arr) noexcept {
    using Result = typename LazyArray<Func, 1>::Result;
    if (arr.sizes[0] == 0)
        return Result{};
    Result res = arr(0);
    for (ssize_t ind = 1; ind < arr.sizes[0]; ind++)
        res = std::max(res, arr(ind));
    return res;
}

/// Sum all.
template <typename Func, size_t Rank>
[[gnu::always_inline]] constexpr auto sum(
    const LazyArray<Func, Rank>& arr) noexcept {
    using Result = typename LazyArray<Func, Rank>::Result;
    Result res = {};
    for (auto ind : ArrayRange(arr.sizes))
        res += arr(ind);
    return res;
}

/// Vector dot.
template <typename Func0, typename Func1>
[[gnu::always_inline]] constexpr auto dot(
    const LazyArray<Func0, 1>& arr0,
    const LazyArray<Func1, 1>& arr1) noexcept {
    ASSERT(arr0.sizes == arr1.sizes);
    using Result0 = typename LazyArray<Func0, 1>::Result;
    using Result1 = typename LazyArray<Func1, 1>::Result;
    using Result = decltype(Result0() * Result1());
    Result res = {};
    for (ssize_t ind = 0; ind < arr0.sizes[0]; ind++)
        res += arr0(ind) * arr1(ind);
    return res;
}

/// Matrix dot.
template <typename Func0, typename Func1>
[[gnu::always_inline]] constexpr auto dot(
    const LazyArray<Func0, 2>& arr0,
    const LazyArray<Func1, 2>& arr1) noexcept {
    ASSERT(arr0.sizes[1] == arr1.sizes[0]);
    using Result0 = typename LazyArray<Func0, 2>::Result;
    using Result1 = typename LazyArray<Func1, 2>::Result;
    using Result = decltype(Result0() * Result1());
    return LazyArray(
        [&](auto ind) constexpr noexcept {
            Result res = {};
            for (ssize_t inner = 0; inner < arr0.sizes[1]; inner++)
                res += arr0(ind[0], inner) * arr1(inner, ind[1]);
            return res;
        },
        ArrayIndex{
            arr0.sizes[0], //
            arr1.sizes[1]});
}

/// Generalized tensor dot.
template <
    typename Func0,
    typename Func1,
    size_t Rank0,
    size_t Rank1,
    size_t OmitRank>
[[gnu::always_inline]] constexpr auto dot(
    const LazyArray<Func0, Rank0>& arr0,
    const LazyArray<Func1, Rank1>& arr1,
    const ArrayIndex<OmitRank>& dims0,
    const ArrayIndex<OmitRank>& dims1) noexcept {
    using Result0 = typename LazyArray<Func0, Rank0>::Result;
    using Result1 = typename LazyArray<Func1, Rank1>::Result;
    using Result = decltype(Result0() * Result1());
    auto [kept0, omitted0] = arr0.sizes.omit(dims0);
    auto [kept1, omitted1] = arr1.sizes.omit(dims1);
    ASSERT(omitted0 == omitted1);
    ArrayIndex<OmitRank> omitted = omitted0;
    return LazyArray(
        [&, omitted ](auto ind) constexpr noexcept {
            Result res = {};
            auto itr = ind.begin();
            auto ind0 = ArrayIndex<Rank0>(itr, dims0);
            auto ind1 = ArrayIndex<Rank1>(itr, dims1);
            for (auto inner : ArrayRange(omitted)) {
                for (size_t dim = 0; dim < OmitRank; dim++) {
                    ind0[dims0[dim]] = inner[dim];
                    ind1[dims1[dim]] = inner[dim];
                }
                res += arr0(ind0) * arr1(ind1);
            }
            return res;
        },
        kept0.join(kept1));
}

/// Outer product.
template <typename Func0, typename Func1, size_t Rank0, size_t Rank1>
[[gnu::always_inline]] constexpr auto outer(
    const LazyArray<Func0, Rank0>& arr0,
    const LazyArray<Func1, Rank1>& arr1) noexcept {
    return LazyArray(
        [&](auto ind) constexpr noexcept {
            return arr0(ind.template head<Rank0>()) *
                   arr1(ind.template tail<Rank1>());
        },
        arr0.sizes.join(arr1.sizes));
}

/// Trace.
template <typename Func>
[[gnu::always_inline]] constexpr auto trace(
    const LazyArray<Func, 2>& arr) noexcept {
    using Result = typename LazyArray<Func, 2>::Result;
    Result res = 0;
    for (ssize_t ind = 0; ind < std::min(arr.sizes[0], arr.sizes[1]); ind++)
        res += arr(ind, ind);
    return res;
}

/// Transpose dimensions.
template <typename Func>
[[gnu::always_inline]] constexpr auto transpose(
    const LazyArray<Func, 2>& arr) noexcept {
    return LazyArray(
        [&](auto ind) constexpr noexcept {
            return arr(ArrayIndex{ind[1], ind[0]});
        },
        ArrayIndex{
            arr.sizes[1], //
            arr.sizes[0]});
}

/// Levi-Civita tensor.
template <size_t Rank>
[[gnu::always_inline]] constexpr auto levi_civita() noexcept {
    ArrayIndex<Rank> sizes;
    std::fill(sizes.begin(), sizes.end(), Rank);
    return LazyArray(
        [](auto ind) constexpr noexcept { return ind.levi_civita(); }, sizes);
}

/// Levi-Civita tensor.
template <size_t Rank>
[[gnu::always_inline]] constexpr auto levi_civita(
    const ArrayIndex<Rank>& sizes) noexcept {
    return LazyArray(
        [](auto ind) constexpr noexcept { return ind.levi_civita(); }, sizes);
}

#if 0
/// Reorder dimensions.
template <typename Func, size_t Rank>
[[gnu::always_inline]] constexpr auto reorder(
    const LazyArray<Func, Rank>& arr,
    const std::type_identity_t<ArrayIndex<Rank>>& ord) noexcept
    requires(Rank > 2) {
    ArrayIndex<Rank> sizes;
    for (size_t dim = 0; dim < Rank; dim++)
        sizes[dim] = arr.sizes[ord[dim]];
    return LazyArray(
        [&](auto ind) constexpr noexcept {
            ArrayIndex<Rank> reorder;
            for (size_t dim = 0; dim < Rank; dim++)
                reorder[dim] = ind[ord[dim]];
            return arr(reorder);
        },
        sizes);
}
#endif

} // namespace pre

#include "Lazy_operators.inl"
