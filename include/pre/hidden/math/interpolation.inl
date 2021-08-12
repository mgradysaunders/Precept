/*-*- C++ -*-*/
#pragma once

namespace pre {

/// Linear interpolation.
///
/// \param[in] t   Factor.
/// \param[in] p0  Control point at \f$ t = 0 \f$.
/// \param[in] p1  Control point at \f$ t = 1 \f$.
///
template <std::floating_point Float, typename Control>
constexpr auto lerp(Float t, const Control& p0, const Control& p1) noexcept {
    return (1 - t) * p0 + t * p1;
}

/// Ease in and out.
template <std::floating_point Float>
[[gnu::always_inline]] constexpr Float ease(Float t) noexcept {
    return t * t * (3 - 2 * t);
}

/// Ease in.
template <std::floating_point Float>
[[gnu::always_inline]] constexpr Float ease_in(Float t) noexcept {
    return t * t;
}

/// Ease out.
template <std::floating_point Float>
[[gnu::always_inline]] constexpr Float ease_out(Float t) noexcept {
    return t * (2 - t);
}

/// Elastic ease in.
template <std::floating_point Float>
[[gnu::always_inline]] inline Float elastic_in(Float t) noexcept {
    if (0 < t and t < 1) {
        constexpr Float a = M_PI * 2 / 3 * 8;
        constexpr Float b = M_PI_2;
        Float u = 8 * t - 8;
        return -Float(t * std::exp2(u) * std::sin(a * u - b));
    }
    else
        return t;
}

/// Elastic ease out.
template <std::floating_point Float>
[[gnu::always_inline]] inline Float elastic_out(Float t) noexcept {
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
template <std::floating_point Float>
constexpr Float smoothstep(Float x, Float x0, Float x1) noexcept {
    Float t = (x - x0) / (x1 - x0);
    t = pre::max(t, Float(0));
    t = pre::min(t, Float(1));
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
template <std::floating_point Float, typename Control>
constexpr auto hermite(
    Float t,
    const Control& p0,
    const Control& m0,
    const Control& m1,
    const Control& p1) noexcept {
    Float s = t - 1;
    Float h00 = s * s * (1 + 2 * t), h10 = s * s * t;
    Float h01 = t * t * (3 - 2 * t), h11 = t * t * s;
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
template <std::floating_point Float, typename Control>
constexpr auto hermite_deriv(
    Float t,
    const Control& p0,
    const Control& m0,
    const Control& m1,
    const Control& p1) noexcept {
    Float g00 = 6 * t * (t - 1);
    Float g10 = 3 * t * t - 4 * t + 1;
    Float g11 = 3 * t * t - 2 * t;
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
template <std::floating_point Float, typename Control>
constexpr auto catmull_rom(
    Float t,
    const Control& pprev,
    const Control& p0,
    const Control& p1,
    const Control& pnext) noexcept {
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
template <std::floating_point Float, typename Control>
constexpr auto catmull_rom_deriv(
    Float t,
    const Control& pprev,
    const Control& p0,
    const Control& p1,
    const Control& pnext) noexcept {
    return hermite_deriv(t, p0, p1 - pprev, pnext - p0, p1);
}

template <std::floating_point Float, typename Control>
constexpr auto bezier(
    Float t, const Control& a, const Control& b, const Control& c) noexcept {
    Float u = 1 - t;
    return u * u * a + Float(2) * u * t * b + t * t * c;
}

template <std::floating_point Float, typename Control>
constexpr auto bezier(
    Float t,
    const Control& a,
    const Control& b,
    const Control& c,
    const Control& d) noexcept {
    Float u = 1 - t;
    Float u2 = u * u, u3 = u2 * u;
    Float t2 = t * t, t3 = t2 * t;
    return u3 * a + Float(3) * (u2 * t * b + u * t2 * c) + t3 * d;
}

} // namespace pre
