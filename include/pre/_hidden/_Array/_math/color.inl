/*-*- C++ -*-*/
#pragma once

namespace pre {

/// Encode linear RGB as sRGB.
///
/// \par Expression
/// \f[
///      \begin{bmatrix} r' \\ g' \\ b' \end{bmatrix} \gets
///      \begin{bmatrix}
///          \operatorname{srgbenc}(r)
///      \\  \operatorname{srgbenc}(g)
///      \\  \operatorname{srgbenc}(b)
///      \end{bmatrix}
/// \f]
///
template <std::floating_point Float>
inline Array<Float, 3> srgbenc(const Array<Float, 3>& v) noexcept {
    return {srgbenc(v[0]), srgbenc(v[1]), srgbenc(v[2])};
}

/// Encode linear RGB+A as sRGB+A.
///
/// \par Expression
/// \f[
///      \begin{bmatrix} r' \\ g' \\ b' \\ \alpha' \end{bmatrix} \gets
///      \begin{bmatrix}
///          \operatorname{srgbenc}(r)
///      \\  \operatorname{srgbenc}(g)
///      \\  \operatorname{srgbenc}(b)
///      \\  \alpha
///      \end{bmatrix}
/// \f]
///
template <std::floating_point Float>
inline Array<Float, 4> srgbenc(const Array<Float, 4>& v) noexcept {
    return {srgbenc(v[0]), srgbenc(v[1]), srgbenc(v[2]), v[3]};
}

/// Decode linear RGB from sRGB.
///
/// \par Expression
/// \f[
///      \begin{bmatrix} r' \\ g' \\ b' \end{bmatrix} \gets
///      \begin{bmatrix}
///          \operatorname{srgbdec}(r)
///      \\  \operatorname{srgbdec}(g)
///      \\  \operatorname{srgbdec}(b)
///      \end{bmatrix}
/// \f]
///
template <std::floating_point Float>
inline Array<Float, 3> srgbdec(const Array<Float, 3>& v) noexcept {
    return {srgbdec(v[0]), srgbdec(v[1]), srgbdec(v[2])};
}

/// Decode linear RGB+A from sRGB+A.
///
/// \par Expression
/// \f[
///      \begin{bmatrix} r' \\ g' \\ b' \\ \alpha' \end{bmatrix} \gets
///      \begin{bmatrix}
///          \operatorname{srgbdec}(r)
///      \\  \operatorname{srgbdec}(g)
///      \\  \operatorname{srgbdec}(b)
///      \\  \alpha
///      \end{bmatrix}
/// \f]
///
template <std::floating_point Float>
inline Array<Float, 4> srgbdec(const Array<Float, 4>& v) noexcept {
    return {srgbdec(v[0]), srgbdec(v[1]), srgbdec(v[2]), v[3]};
}

/// Encode linear RGB as sRGB with Hable tonemapping.
///
/// \par Expression
/// \f[
///      \begin{bmatrix} r' \\ g' \\ b' \end{bmatrix} \gets
///      \begin{bmatrix}
///          \operatorname{srgbenc}_{\text{Hable}}(r)
///      \\  \operatorname{srgbenc}_{\text{Hable}}(g)
///      \\  \operatorname{srgbenc}_{\text{Hable}}(b)
///      \end{bmatrix}
/// \f]
///
/// \see
/// [This article][1] by John Hable.
/// [1]: http://filmicworlds.com/blog/filmic-tonemapping-operators/
///
template <std::floating_point Float>
inline Array<Float, 3> srgbenc_hable(const Array<Float, 3>& v) noexcept {
    return {srgbenc_hable(v[0]), srgbenc_hable(v[1]), srgbenc_hable(v[2])};
}

/// Encode linear RGB+A as sRGB+A with Hable tonemapping.
///
/// \par Expression
/// \f[
///      \begin{bmatrix} r' \\ g' \\ b' \\ \alpha' \end{bmatrix} \gets
///      \begin{bmatrix}
///          \operatorname{srgbenc}_{\text{Hable}}(r)
///      \\  \operatorname{srgbenc}_{\text{Hable}}(g)
///      \\  \operatorname{srgbenc}_{\text{Hable}}(b)
///      \\  \alpha
///      \end{bmatrix}
/// \f]
///
/// \see
/// [This article][1] by John Hable.
/// [1]: http://filmicworlds.com/blog/filmic-tonemapping-operators/
///
template <std::floating_point Float>
inline Array<Float, 4> srgbenc_hable(const Array<Float, 4>& v) noexcept {
    return {srgbenc_hable(v[0]), srgbenc_hable(v[1]), srgbenc_hable(v[2]),
            v[3]};
}

/// Encode linear RGB as sRGB with Hejl/Burgess-Dawson tonemapping.
///
/// \par Expression
/// \f[
///      \begin{bmatrix} r' \\ g' \\ b' \end{bmatrix} \gets
///      \begin{bmatrix}
///          \operatorname{srgbenc}_{\text{Hejl/Burgess}}(r)
///      \\  \operatorname{srgbenc}_{\text{Hejl/Burgess}}(g)
///      \\  \operatorname{srgbenc}_{\text{Hejl/Burgess}}(b)
///      \end{bmatrix}
/// \f]
///
/// \see
/// [This article][1] by John Hable.
/// [1]: http://filmicworlds.com/blog/filmic-tonemapping-operators/
///
template <std::floating_point Float>
inline Array<Float, 3> srgbenc_hejl_burgess(
        const Array<Float, 3>& v) noexcept {
    return {srgbenc_hejl_burgess(v[0]), srgbenc_hejl_burgess(v[1]),
            srgbenc_hejl_burgess(v[2])};
}

/// Encode linear RGB+A as sRGB+A with Hejl/Burgess-Dawson tonemapping.
///
/// \par Expression
/// \f[
///      \begin{bmatrix} r' \\ g' \\ b' \\ \alpha' \end{bmatrix} \gets
///      \begin{bmatrix}
///          \operatorname{srgbenc}_{\text{Hejl/Burgess}}(r)
///      \\  \operatorname{srgbenc}_{\text{Hejl/Burgess}}(g)
///      \\  \operatorname{srgbenc}_{\text{Hejl/Burgess}}(b)
///      \\  \alpha
///      \end{bmatrix}
/// \f]
///
/// \see
/// [This article][1] by John Hable.
/// [1]: http://filmicworlds.com/blog/filmic-tonemapping-operators/
///
template <std::floating_point Float>
inline Array<Float, 4> srgbenc_hejl_burgess(
        const Array<Float, 4>& v) noexcept {
    return {srgbenc_hejl_burgess(v[0]), srgbenc_hejl_burgess(v[1]),
            srgbenc_hejl_burgess(v[2]), v[3]};
}

/// XYZ triple to RGB triple.
///
/// \note
/// This uses standard CIE parameters:
/// - \f$ C_r = (0.7350, 0.2650) \f$,
/// - \f$ C_g = (0.2740, 0.7170) \f$,
/// - \f$ C_b = (0.1670, 0.0090) \f$, and
/// - \f$ W = (1, 1, 1) \f$.
///
/// \see
/// [Bruce Lindbloom's page][1].
/// [1]: http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
///
template <std::floating_point Float>
inline Array<Float, 3> xyz_to_rgb(const Array<Float, 3>& v) noexcept {
    Array<Float, 3, 3> m = {
            Float(+2.3706743), Float(-0.9000405), Float(-0.4706338), //
            Float(-0.5138850), Float(+1.4253036), Float(+0.0885814), //
            Float(+0.0052982), Float(-0.0146949), Float(+1.0093968)};
    return dot(m, v);
}

/// RGB triple to XYZ triple.
///
/// \note
/// This uses standard CIE parameters:
/// - \f$ C_r = (0.7350, 0.2650) \f$,
/// - \f$ C_g = (0.2740, 0.7170) \f$,
/// - \f$ C_b = (0.1670, 0.0090) \f$, and
/// - \f$ W = (1, 1, 1) \f$.
///
/// \see
/// [Bruce Lindbloom's page][1].
/// [1]: http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
///
template <std::floating_point Float>
inline Array<Float, 3> rgb_to_xyz(const Array<Float, 3>& v) noexcept {
    Array<Float, 3, 3> m = {
            Float(0.4887180), Float(0.3106803), Float(0.2006017), //
            Float(0.1762044), Float(0.8129847), Float(0.0108109), //
            Float(0.0000000), Float(0.0102048), Float(0.9897952)};
    return dot(m, v);
}

/// RGB to XYZ conversion matrix.
///
/// \param[in] cr
/// xy pair of reference red.
///
/// \param[in] cg
/// xy pair of reference green.
///
/// \param[in] cb
/// xy pair of reference blue.
///
/// \param[in] w
/// XYZ triple of reference white.
///
/// \see
/// [Bruce Lindbloom's page][1].
/// [1]: http://brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
///
template <std::floating_point Float>
inline Array<Float, 3, 3> rgb_to_xyz(
        const Array<Float, 2>& cr,
        const Array<Float, 2>& cg,
        const Array<Float, 2>& cb,
        const Array<Float, 3>& w) noexcept {
    // Conversion matrix.
    Array<Float, 3, 3> m;

    // Temp matrix.
    Array<Float, 3, 3> a = {
            cr[0] / cr[1],
            cg[0] / cg[1],
            cb[0] / cb[1], //
            Float(1),
            Float(1),
            Float(1), //
            (Float(1) - cr[0] - cr[1]) / cr[1],
            (Float(1) - cg[0] - cg[1]) / cg[1],
            (Float(1) - cb[0] - cb[1]) / cb[1]};

    // Temp matrix cofactors.
    Array<Float, 3, 3> acof;
    acof[0] = cross(a[1], a[2]);
    acof[1] = cross(a[2], a[0]);
    acof[2] = cross(a[0], a[1]);

    // Temp matrix inverse times reference white.
    Array<Float, 3> s = dot(w, acof) / acof(1).sum();

    // Initialize conversion matrix.
    m[0] = a[0] * s;
    m[1] = a[1] * s;
    m[2] = a[2] * s;
    return m;
}

/// XYZ triple to xyY triple.
///
/// \par Expression
/// \f[
///      \begin{bmatrix} x \\ y \\ Y \end{bmatrix} \gets
///      \begin{bmatrix}
///          X / (X + Y + Z)
///      \\  Y / (X + Y + Z)
///      \\  Y
///      \end{bmatrix}
/// \f]
///
template <std::floating_point Float>
inline Array<Float, 3> xyz_to_xyy(const Array<Float, 3>& v) noexcept {
    return {v[0] / v.sum(), v[1] / v.sum(), v[1]};
}

/// xyY triple to XYZ triple.
///
/// \par Expression
/// \f[
///      \begin{bmatrix} X \\ Y \\ Z \end{bmatrix} \gets
///      Y
///      \begin{bmatrix}
///          x / y
///      \\  1
///      \\  (1 - x - y) / y
///      \end{bmatrix}
/// \f]
///
template <std::floating_point Float>
inline Array<Float, 3> xyy_to_xyz(const Array<Float, 3>& v) noexcept {
    return {v[2] / v[1] * v[0], v[2], v[2] / v[1] * (1 - v[0] - v[1])};
}

/// XYZ triple to Lab triple.
///
/// \par Expression
/// \f[
///      \begin{bmatrix} L \\ a \\ b \end{bmatrix} \gets
///      \begin{bmatrix}
///          116 f(Y) - 16
///      \\  500 (f(X) - f(Y))
///      \\  200 (f(Y) - f(Z))
///      \end{bmatrix}
/// \f]
/// where
/// \f[
///      f(t) =
///      \begin{cases}
///          t^{1/3}                     & t >   216 / 24389
///      \\  ((24389 / 27) t + 16) / 116 & t \le 216 / 24389
///      \end{cases}
/// \f]
///
/// \see
/// [Bruce Lindbloom's page][1].
/// [1]: http://brucelindbloom.com/index.html?Eqn_XYZ_to_Lab.html
///
template <std::floating_point Float>
inline Array<Float, 3> xyz_to_lab(const Array<Float, 3>& v) noexcept {
    auto f = [](Float t) {
        if (t > Float(216) / Float(24389))
            return pre::cbrt(t);
        else
            return (t * (Float(24389) / Float(27)) + Float(16)) / Float(116);
    };
    return {Float(116) * f(v[1]) - Float(16), //
            Float(500) * (f(v[0]) - f(v[1])), //
            Float(200) * (f(v[1]) - f(v[2]))};
}

/// Lab triple to XYZ triple.
///
/// \par Expression
/// \f[
///      \begin{bmatrix} X \\ Y \\ Z \end{bmatrix} \gets
///      \begin{bmatrix}
///          f^{-1}(f_X)
///      \\  f^{-1}(f_Y)
///      \\  f^{-1}(f_Z)
///      \end{bmatrix}
/// \f]
/// where
/// \f[
///      \begin{aligned}
///          f_Y &= (L + 16) / 116
///      \\  f_X &= f_Y + a / 500
///      \\  f_Z &= f_Y - b / 200
///      \end{aligned}
/// \f]
/// and where
/// \f[
///      f^{-1}(t) =
///      \begin{cases}
///          t^3                       & t^3 >   216 / 24389
///      \\  (116 t - 16) (27 / 24389) & t^3 \le 216 / 24389
///      \end{cases}
/// \f]
///
/// \see
/// [Bruce Lindbloom's page][1].
/// [1]: http://brucelindbloom.com/index.html?Eqn_Lab_to_XYZ.html
///
template <std::floating_point Float>
inline Array<Float, 3> lab_to_xyz(const Array<Float, 3>& v) noexcept {
    auto finv = [](Float t) {
        if (nthpow(t, 3) > Float(216) / Float(24389))
            return nthpow(t, 3);
        else
            return (t * Float(116) - Float(16)) * (Float(27) / Float(24389));
    };
    Float fy = (v[0] + Float(16)) / Float(116);
    Float fx = fy + v[1] / Float(500);
    Float fz = fy - v[2] / Float(200);
    return {finv(fx), finv(fy), finv(fz)};
}

/// XYZ triple to Luv triple.
///
/// \param[in] w  XYZ triple of reference white.
/// \param[in] v  XYZ triple.
///
template <std::floating_point Float>
inline Array<Float, 3> xyz_to_luv(
        const Array<Float, 3>& w, const Array<Float, 3>& v) noexcept {
    Float t0 = nthpow(Float(6) / Float(29), 3);
    Float t = v[1] / w[1];
    Float l = t <= t0 ? t / t0 : 116 * pre::cbrt(t) - 16;
    Array<Float, 2> wp = {4 * w[0], 9 * w[1]};
    Array<Float, 2> vp = {4 * v[0], 9 * v[1]};
    wp /= w[0] + 15 * w[1] + 3 * w[2];
    vp /= v[0] + 15 * v[1] + 3 * v[2];
    return {l, 13 * l * (vp[0] - wp[0]), 13 * l * (vp[1] - wp[1])};
}

/// Luv triple to XYZ triple.
///
/// \param[in] w  XYZ triple of reference white.
/// \param[in] v  Luv triple.
///
template <std::floating_point Float>
inline Array<Float, 3> luv_to_xyz(
        const Array<Float, 3>& w, const Array<Float, 3>& v) noexcept {
    Float y = v[0] <= 8 ? v[0] * nthpow(Float(6) / Float(29), 3)
                        : nthpow((v[0] + 16) / 116, 3);
    if (!(y > 0))
        return {};
    y *= w[1];
    Array<Float, 2> wp = {
            4 * w[0] / (w[0] + 15 * w[1] + 3 * w[2]),
            9 * w[1] / (w[0] + 15 * w[1] + 3 * w[2])};
    Array<Float, 2> vp = {
            wp[0] + v[1] / (13 * v[0]), wp[1] + v[2] / (13 * v[0])};
    if (!pre::isfinite(vp).all())
        return {};
    else
        return {y / (4 * vp[1]) * (9 * vp[0]), y,
                y / (4 * vp[1]) * (12 - 3 * vp[0] - 20 * vp[1])};
}

/// RGB triple to Lab triple.
template <std::floating_point Float>
inline Array<Float, 3> rgb_to_lab(const Array<Float, 3>& v) noexcept {
    return xyz_to_lab(rgb_to_xyz(v));
}

/// Lab triple to RGB triple.
template <std::floating_point Float>
inline Array<Float, 3> lab_to_rgb(const Array<Float, 3>& v) noexcept {
    return xyz_to_rgb(lab_to_xyz(v));
}

/// Mix RGB.
///
/// \param[in] t   Factor in \f$ [0, 1] \f$.
/// \param[in] v0  RGB triple at \f$ t = 0 \f$.
/// \param[in] v1  RGB triple at \f$ t = 1 \f$.
///
template <std::floating_point Float>
inline Array<Float, 3> mix_rgb(
        Float t,
        const Array<Float, 3>& v0,
        const Array<Float, 3>& v1) noexcept {
    return lab_to_rgb(lerp(t, rgb_to_lab(v0), rgb_to_lab(v1)));
}

/// Mix RGB+A.
///
/// \param[in] t   Factor in \f$ [0, 1] \f$.
/// \param[in] v0  RGB+A quadruple at \f$ t = 0 \f$.
/// \param[in] v1  RGB+A quadruple at \f$ t = 1 \f$.
///
template <std::floating_point Float>
inline Array<Float, 4> mix_rgb(
        Float t,
        const Array<Float, 4>& v0,
        const Array<Float, 4>& v1) noexcept {
    Array<Float, 4> v = mix_rgb(t, Array<Float, 3>(v0), Array<Float, 3>(v1));
    v[3] = lerp(t, v0[3], v1[3]);
    return v;
}

// TODO color literals

} // namespace pre
