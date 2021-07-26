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
template <std::floating_point Float>
inline Float srgbenc(Float v) {
    if (v <= Float(0.0031308)) {
        return Float(12.92) * v;
    }
    else {
        return Float(1.055) * pre::pow(v, Float(1) / Float(2.4)) -
               Float(0.055);
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
template <std::floating_point Float>
inline Float srgbdec(Float v) {
    if (v <= Float(0.04045)) {
        return v / Float(12.92);
    }
    else {
        return pre::pow((v + Float(0.055)) / Float(1.055), Float(2.4));
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
template <std::floating_point Float>
inline Float srgbenc_hable(Float v) {
    return srgbenc(
            (v * (v * Float(5.50710) + Float(0.91785)) + Float(0.036714)) /
                    (v * (v * Float(3.99336) + Float(6.65560)) +
                     Float(0.399336)) -
            Float(0.6119) / Float(6.6556));
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
template <std::floating_point Float>
inline Float srgbenc_hejl_burgess(Float v) {
    v = pre::max(v - Float(0.004), Float(0));
    return (v * (v * Float(6.2) + Float(0.5))) /
           (v * (v * Float(6.2) + Float(1.7)) + Float(0.06));
}

/// Fit of CIE 1931 X by Wyman et al.
///
/// \param[in] lambda
/// Wavelength in micrometers.
///
/// \see [This publication][1] by Wyman, Sloan, and Shirley.
/// [1]: http://jcgt.org/published/0002/02/01/
///
template <std::floating_point Float>
inline Float wymanx(Float lambda) {
    Float t1 = lambda - Float(0.4420);
    Float t2 = lambda - Float(0.5998);
    Float t3 = lambda - Float(0.5011);
    t1 *= pre::signbit(t1) ? Float(62.4) : Float(37.4);
    t2 *= pre::signbit(t2) ? Float(26.4) : Float(32.3);
    t3 *= pre::signbit(t3) ? Float(49.0) : Float(38.2);
    return Float(0.362) * pre::exp(Float(-0.5) * t1 * t1) +
           Float(1.056) * pre::exp(Float(-0.5) * t2 * t2) -
           Float(0.065) * pre::exp(Float(-0.5) * t3 * t3);
}

/// Fit of CIE 1931 Y by Wyman et al.
///
/// \param[in] lambda
/// Wavelength in micrometers.
///
/// \see [This publication][1] by Wyman, Sloan, and Shirley.
/// [1]: http://jcgt.org/published/0002/02/01/
///
template <std::floating_point Float>
inline Float wymany(Float lambda) {
    Float t1 = lambda - Float(0.5688);
    Float t2 = lambda - Float(0.5309);
    t1 *= pre::signbit(t1) ? Float(21.3) : Float(24.7);
    t2 *= pre::signbit(t2) ? Float(61.3) : Float(32.2);
    return Float(0.821) * pre::exp(Float(-0.5) * t1 * t1) +
           Float(0.286) * pre::exp(Float(-0.5) * t2 * t2);
}

/// Fit of CIE 1931 Z by Wyman et al.
///
/// \param[in] lambda
/// Wavelength in micrometers.
///
/// \see [This publication][1] by Wyman, Sloan, and Shirley.
/// [1]: http://jcgt.org/published/0002/02/01/
///
template <std::floating_point Float>
inline Float wymanz(Float lambda) {
    Float t1 = lambda - Float(0.4370);
    Float t2 = lambda - Float(0.4590);
    t1 *= pre::signbit(t1) ? Float(84.5) : Float(27.8);
    t2 *= pre::signbit(t2) ? Float(38.5) : Float(72.5);
    return Float(1.217) * pre::exp(Float(-0.5) * t1 * t1) +
           Float(0.681) * pre::exp(Float(-0.5) * t2 * t2);
}

/// Planck's law.
///
/// Planck's law of blackbody radiation is
/// \f[
///     b(Float,\lambda) =
///         \frac{1}{\lambda^5}
///         \frac{2hc^2}{e^{\frac{hc}{kT\lambda}}-1}
/// \f]
/// where, by typical conventions, \f$ Float \f$ is temperature
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
template <std::floating_point Float>
inline Float planck(Float t, Float lambda) {
    if (!(lambda > 0)) {
        return 0;
    }
    else {
        constexpr Float c0 = Float(1.19104290768681554502861912e+02L);
        constexpr Float c1 = Float(1.43877729954300303744214349e+04L);
        Float lambda2 = lambda * lambda;
        Float lambda4 = lambda2 * lambda2;
        Float lambda5 = lambda4 * lambda;
        return c0 / (lambda5 * pre::expm1(c1 / (t * lambda)));
    }
}

} // namespace pre
