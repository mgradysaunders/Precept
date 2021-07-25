/*-*- C++ -*-*/
#pragma once

namespace pre {

/// Encode linear RGB as sRGB.
///
/// \par Expression
/// \f[
///     \text{srgbenc}(v) =
///     \begin{cases}
///         12.92 v                 & v \le 0.0031308
///     \\  1.055 v^{1/2.4} - 0.055 & v >   0.0031308
///     \end{cases}
/// \f]
///
template <std::floating_point T>
inline T srgbenc(T v) {
    if (v <= T(0.0031308)) {
        return T(12.92) * v;
    }
    else {
        return T(1.055) * pre::pow(v, T(1) / T(2.4)) - T(0.055);
    }
}

/// Decode linear RGB from sRGB.
///
/// \par Expression
/// \f[
///     \text{srgbdec}(v) =
///     \begin{cases}
///         v / 12.92                   & v \le 0.04045
///     \\  ((v + 0.055) / 1.055)^{2.4} & v >   0.04045
///     \end{cases}
/// \f]
///
template <std::floating_point T>
inline T srgbdec(T v) {
    if (v <= T(0.04045)) {
        return v / T(12.92);
    }
    else {
        return pre::pow((v + T(0.055)) / T(1.055), T(2.4));
    }
}

/// Encode linear RGB as sRGB with Hable tonemapping.
///
/// \par Expression
/// \f[
///     \text{srgbenc}_{\text{Hable}}(v) =
///     \text{srgbenc}\left(
///     \frac{5.50710 v^2 + 0.91785 v + 0.036714}
///          {3.99336 v^2 + 6.65560 v + 0.399336} -
///     \frac{0.61190}{6.65560}\right)
/// \f]
///
/// \see
/// [This article][1] by John Hable.
/// [1]: http://filmicworlds.com/blog/filmic-tonemapping-operators/
///
template <std::floating_point T>
inline T srgbenc_hable(T v) {
    return srgbenc(
            (v * (v * T(5.50710) + T(0.91785)) + T(0.036714)) /
                    (v * (v * T(3.99336) + T(6.65560)) + T(0.399336)) -
            T(0.6119) / T(6.6556));
}

/// Encode linear RGB as sRGB with Hejl/Burgess-Dawson tonemapping.
///
/// \par Expression
/// \f[
///     \text{srgbenc}_{\text{Hejl/Burgess}}(v) =
///     \frac{6.2 \max(v - 0.004, 0)^2 + 0.5 \max(v - 0.004, 0)}
///          {6.2 \max(v - 0.004, 0)^2 + 1.7 \max(v - 0.004, 0) + 0.06}
/// \f]
///
/// \see
/// [This article][1] by John Hable.
/// [1]: http://filmicworlds.com/blog/filmic-tonemapping-operators/
///
template <std::floating_point T>
inline T srgbenc_hejl_burgess(T v) {
    v = pre::max(v - T(0.004), T(0));
    return (v * (v * T(6.2) + T(0.5))) / (v * (v * T(6.2) + T(1.7)) + T(0.06));
}

/// Fit of CIE 1931 X by Wyman et al.
///
/// \param[in] lambda
/// Wavelength in micrometers.
///
/// \see [This publication][1] by Wyman, Sloan, and Shirley.
/// [1]: http://jcgt.org/published/0002/02/01/
///
template <std::floating_point T>
inline T wymanx(T lambda) {
    T t1 = lambda - T(0.4420);
    T t2 = lambda - T(0.5998);
    T t3 = lambda - T(0.5011);
    t1 *= pre::signbit(t1) ? T(62.4) : T(37.4);
    t2 *= pre::signbit(t2) ? T(26.4) : T(32.3);
    t3 *= pre::signbit(t3) ? T(49.0) : T(38.2);
    return T(0.362) * pre::exp(T(-0.5) * t1 * t1) +
           T(1.056) * pre::exp(T(-0.5) * t2 * t2) -
           T(0.065) * pre::exp(T(-0.5) * t3 * t3);
}

/// Fit of CIE 1931 Y by Wyman et al.
///
/// \param[in] lambda
/// Wavelength in micrometers.
///
/// \see [This publication][1] by Wyman, Sloan, and Shirley.
/// [1]: http://jcgt.org/published/0002/02/01/
///
template <std::floating_point T>
inline T wymany(T lambda) {
    T t1 = lambda - T(0.5688);
    T t2 = lambda - T(0.5309);
    t1 *= pre::signbit(t1) ? T(21.3) : T(24.7);
    t2 *= pre::signbit(t2) ? T(61.3) : T(32.2);
    return T(0.821) * pre::exp(T(-0.5) * t1 * t1) +
           T(0.286) * pre::exp(T(-0.5) * t2 * t2);
}

/// Fit of CIE 1931 Z by Wyman et al.
///
/// \param[in] lambda
/// Wavelength in micrometers.
///
/// \see [This publication][1] by Wyman, Sloan, and Shirley.
/// [1]: http://jcgt.org/published/0002/02/01/
///
template <std::floating_point T>
inline T wymanz(T lambda) {
    T t1 = lambda - T(0.4370);
    T t2 = lambda - T(0.4590);
    t1 *= pre::signbit(t1) ? T(84.5) : T(27.8);
    t2 *= pre::signbit(t2) ? T(38.5) : T(72.5);
    return T(1.217) * pre::exp(T(-0.5) * t1 * t1) +
           T(0.681) * pre::exp(T(-0.5) * t2 * t2);
}

/// Planck's law.
///
/// Planck's law of blackbody radiation is
/// \f[
///     b(T,\lambda) =
///         \frac{1}{\lambda^5}
///         \frac{2hc^2}{e^{\frac{hc}{kT\lambda}}-1}
/// \f]
/// where, by typical conventions, \f$ T \f$ is temperature
/// in degrees kelvin and \f$ \lambda \f$ is wavelength in meters. The
/// implementation here takes \f$ \lambda \f$ in _micrometers_ instead of
/// meters, but in the interest of avoiding astronomic values, the output
/// units are \f$ \mathrm{MW}/\mathrm{sr}/\mathrm{m}^{2}/\mu\mathrm{m} \f$.
///
/// \param[in] t
/// Blackbody temperature in degrees kelvin.
///
/// \param[in] lambda
/// Wavelength in micrometers.
///
/// \see Wikipedia's article for [Planck's Law][1].
/// [1]: https://en.wikipedia.org/wiki/Planck%27s_law
///
template <std::floating_point T>
inline T planck(T t, T lambda) {
    if (!(lambda > 0)) {
        return 0;
    }
    else {
        constexpr T c0 = T(1.19104290768681554502861912e+02L);
        constexpr T c1 = T(1.43877729954300303744214349e+04L);
        T lambda2 = lambda * lambda;
        T lambda4 = lambda2 * lambda2;
        T lambda5 = lambda4 * lambda;
        return c0 / (lambda5 * pre::expm1(c1 / (t * lambda)));
    }
}

} // namespace pre
