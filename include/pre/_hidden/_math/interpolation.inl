/*-*- C++ -*-*/
#pragma once

namespace pre {

/// Linear interpolation.
///
/// \param[in] t   Factor.
/// \param[in] p0  Control point at \f$ t = 0 \f$.
/// \param[in] p1  Control point at \f$ t = 1 \f$.
///
template <typename T, typename U>
[[gnu::always_inline]] constexpr auto lerp(T t, const U& p0, const U& p1) {
    return (1 - t) * p0 + t * p1;
}

/// Ease in and out.
template <typename T>
[[gnu::always_inline]] constexpr T ease(T t) {
    return t * t * (3 - 2 * t);
}

/// Ease in.
template <typename T>
[[gnu::always_inline]] constexpr T ease_in(T t) {
    return t * t;
}

/// Ease out.
template <typename T>
[[gnu::always_inline]] constexpr T ease_out(T t) {
    return t * (2 - t);
}

/// Elastic ease in.
template <std::floating_point T>
[[gnu::always_inline]] inline T elastic_in(T t) {
    if (0 < t && t < 1) {
        constexpr T a = M_PI * 2 / 3 * 8;
        constexpr T b = M_PI_2;
        T u = 8 * t - 8;
        return -T(t * std::exp2(u) * std::sin(a * u - b));
    }
    else
        return t;
}

/// Elastic ease out.
template <std::floating_point T>
[[gnu::always_inline]] inline T elastic_out(T t) {
    return 1 - elastic_in(1 - t);
}

/// Smoothstep.
///
/// \param[in] x   Value.
/// \param[in] x0  Value where \f$ t = 0 \f$.
/// \param[in] x1  Value where \f$ t = 1 \f$.
///
/// \par Expression
/// - \f$ t \gets (x - x_0) / (x_1 - x_0) \f$
/// - \f$ t \gets t^2 (3 - 2 t) \f$
///
template <typename T>
[[gnu::always_inline]] constexpr auto smoothstep(T x, T x0, T x1) {
    T t = (x - x0) / (x1 - x0);
    t = pre::max(t, T(0));
    t = pre::min(t, T(1));
    return ease_in_out(t);
}

/// Hermite interpolation.
///
/// \param[in] t   Factor.
/// \param[in] p0  Control point at \f$ t = 0 \f$.
/// \param[in] m0  Slope at \f$ t = 0 \f$.
/// \param[in] m1  Slope at \f$ t = 1 \f$.
/// \param[in] p1  Control point at \f$ t = 1 \f$.
///
/// \see Wikipedia's article for [Cubic Hermite spline][1].
/// [1]: https://en.wikipedia.org/wiki/Cubic_Hermite_spline
///
template <typename T, typename U>
constexpr auto hermite(
        T t, const U& p0, const U& m0, const U& m1, const U& p1) {
    T s = t - 1;
    T h00 = s * s * (1 + 2 * t), h10 = s * s * t;
    T h01 = t * t * (3 - 2 * t), h11 = t * t * s;
    return (h00 * p0 + h10 * m0) + (h01 * p1 + h11 * m1);
}

/// Hermite interpolation deriviatve.
///
/// \param[in] t   Factor.
/// \param[in] p0  Control point at \f$ t = 0 \f$.
/// \param[in] m0  Slope at \f$ t = 0 \f$.
/// \param[in] m1  Slope at \f$ t = 1 \f$.
/// \param[in] p1  Control point at \f$ t = 1 \f$.
///
template <typename T, typename U>
constexpr auto hermite_deriv(
        T t, const U& p0, const U& m0, const U& m1, const U& p1) {
    T g00 = 6 * t * (t - 1);
    T g10 = 3 * t * t - 4 * t + 1;
    T g11 = 3 * t * t - 2 * t;
    return g00 * (p0 - p1) + g10 * m0 + g11 * m1;
}

/// Catmull-Rom interpolation.
///
/// \param[in] t      Factor.
/// \param[in] pprev  Previous control point.
/// \param[in] p0     Control point at \f$ t = 0 \f$.
/// \param[in] p1     Control point at \f$ t = 1 \f$.
/// \param[in] pnext  Next control point.
///
/// \see Wikipedia's article for [Cubic Hermite spline][1].
/// [1]: https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Catmull–Rom_spline
///
template <typename T, typename U>
[[gnu::always_inline]] constexpr auto catmull_rom(
        T t, const U& pprev, const U& p0, const U& p1, const U& pnext) {
    return hermite(t, p0, p1 - pprev, pnext - p0, p1);
}

/// Catmull-Rom interpolation derivative.
///
/// \param[in] t      Factor.
/// \param[in] pprev  Previous control point.
/// \param[in] p0     Control point at \f$ t = 0 \f$.
/// \param[in] p1     Control point at \f$ t = 1 \f$.
/// \param[in] pnext  Next control point.
///
/// \see Wikipedia's article for [Cubic Hermite spline][1].
/// [1]: https://en.wikipedia.org/wiki/Cubic_Hermite_spline#Catmull–Rom_spline
///
template <typename T, typename U>
[[gnu::always_inline]] constexpr auto catmull_rom_deriv(
        T t, const U& pprev, const U& p0, const U& p1, const U& pnext) {
    return hermite_deriv(t, p0, p1 - pprev, pnext - p0, p1);
}

template <typename T, typename U>
[[gnu::always_inline]] constexpr auto bezier(
        T t, const U& a, const U& b, const U& c) {
    T u = T(1) - t;
    return u * u * a + T(2) * u * t * b + t * t * c;
}

template <typename T, typename U>
[[gnu::always_inline]] constexpr auto bezier(
        T t, const U& a, const U& b, const U& c, const U& d) {
    T u = T(1) - t;
    T u2 = u * u, u3 = u2 * u;
    T t2 = t * t, t3 = t2 * t;
    return u3 * a + T(3) * (u2 * t * b + u * t2 * c) + t3 * d;
}

} // namespace pre
