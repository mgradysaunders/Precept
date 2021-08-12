/*-*- C++ -*-*/
#pragma once

namespace pre {

template <typename T, size_t... N>
inline auto min(
        const Array<T, N...>& arrx, const Array<T, N...>& arry) noexcept {
    Array<decltype(pre::min(T(), T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarry = arry.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarry, ++itrres)
        *itrres = pre::min(*itrarrx, *itrarry);
    return res;
}

template <typename T, size_t... N>
inline auto max(
        const Array<T, N...>& arrx, const Array<T, N...>& arry) noexcept {
    Array<decltype(pre::max(T(), T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarry = arry.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarry, ++itrres)
        *itrres = pre::max(*itrarrx, *itrarry);
    return res;
}

template <typename T, size_t... N>
inline auto clamp(
        const Array<T, N...>& arrx,
        const Array<T, N...>& arra,
        const Array<T, N...>& arrb) noexcept {
    Array<decltype(pre::clamp(T(), T(), T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarra = arra.begin();
    auto itrarrb = arrb.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarra, ++itrarrb, ++itrres)
        *itrres = pre::clamp(*itrarrx, *itrarra, *itrarrb);
    return res;
}

template <typename T, size_t... N>
inline auto nthpow(
        const Array<T, N...>& arrx, const Array<int, N...>& arrn) noexcept {
    Array<decltype(pre::nthpow(T(), std::declval<int>())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarrn = arrn.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarrn, ++itrres)
        *itrres = pre::nthpow(*itrarrx, *itrarrn);
    return res;
}

template <typename T, size_t... N>
inline auto rotl(
        const Array<T, N...>& arrval,
        const Array<unsigned, N...>& arrrot) noexcept {
    Array<decltype(pre::rotl(T(), std::declval<unsigned>())), N...> res;
    auto itrarrval = arrval.begin();
    auto itrarrrot = arrrot.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrval, ++itrarrrot, ++itrres)
        *itrres = pre::rotl(*itrarrval, *itrarrrot);
    return res;
}

template <typename T, size_t... N>
inline auto rotr(
        const Array<T, N...>& arrval,
        const Array<unsigned, N...>& arrrot) noexcept {
    Array<decltype(pre::rotr(T(), std::declval<unsigned>())), N...> res;
    auto itrarrval = arrval.begin();
    auto itrarrrot = arrrot.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrval, ++itrarrrot, ++itrres)
        *itrres = pre::rotr(*itrarrval, *itrarrrot);
    return res;
}

template <typename T, size_t... N>
inline auto float_incr(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::float_incr(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::float_incr(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto float_decr(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::float_decr(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::float_decr(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto fast_floor(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::fast_floor(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::fast_floor(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto fast_ceil(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::fast_ceil(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::fast_ceil(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto fast_trunc(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::fast_trunc(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::fast_trunc(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto fast_round(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::fast_round(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::fast_round(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto fract(
        const Array<T, N...>& arrx, Array<int, N...>* arri) noexcept {
    Array<decltype(pre::fract(T(), std::declval<int*>())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrarri = arri->begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrarri, ++itrres)
        *itrres = pre::fract(*itrarrx, &(*itrarri));
    return res;
}

template <typename T, size_t... N>
inline auto sinpi(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::sinpi(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::sinpi(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto cospi(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::cospi(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::cospi(*itrarrx);
    return res;
}

template <typename T, size_t... N>
inline auto erfinv(const Array<T, N...>& arrx) noexcept {
    Array<decltype(pre::erfinv(T())), N...> res;
    auto itrarrx = arrx.begin();
    auto itrres = res.begin();
    for (; itrres < res.end(); ++itrarrx, ++itrres)
        *itrres = pre::erfinv(*itrarrx);
    return res;
}

} // namespace pre
