/*-*- C++ -*-*/
#pragma once

namespace pre {

/// 2-dimensional cross product.
///
/// \par Expression
/// \f[
///      \sum_j \varepsilon_{[i,j]} x_{[j]}
///      \implies
///      \begin{bmatrix}
///          +x_{[1]}
///      \\  -x_{[0]}
///      \end{bmatrix}
/// \f]
///
template <typename Arith>
constexpr Array<Arith, 2> cross(const Array<Arith, 2>& arr) noexcept {
    return {arr[1], -arr[0]};
}

/// 2-dimensional cross product.
///
/// \par Expression
/// \f[
///      \sum_{i,j} \varepsilon_{[i,j]} x_{0[i]} x_{1[j]}
///      \implies
///      x_{0[0]} x_{1[1]} -
///      x_{0[1]} x_{1[0]}
/// \f]
///
template <typename Arith0, typename Arith1>
constexpr auto cross(
    const Array<Arith0, 2>& arr0, const Array<Arith1, 2>& arr1) noexcept {
    return dot(arr0, cross(arr1));
}

/// 3-dimensional cross product.
///
/// \par Expression
/// \f[
///      \sum_k \varepsilon_{[i,j,k]} x_{[k]}
///      \implies
///      \begin{bmatrix}
///           0 & -x_{[2]} & +x_{[1]}
///      \\  +x_{[2]} &  0 & -x_{[0]}
///      \\  -x_{[1]} & +x_{[0]} & 0
///      \end{bmatrix}
/// \f]
///
template <typename Arith>
constexpr Array<Arith, 3, 3> cross(const Array<Arith, 3>& arr) noexcept {
    return {Arith(), -arr[2], +arr[1], //
            +arr[2], Arith(), -arr[0], //
            -arr[1], +arr[0], Arith()};
}

/// 3-dimensional cross product.
///
/// \par Expression
/// \f[
///      \sum_{j,k} \varepsilon_{[i,j,k]} x_{0[j]} x_{1[k]}
///      \implies
///      \begin{bmatrix}
///          x_{0[1]} x_{1[2]} - x_{0[2]} x_{1[1]}
///      \\  x_{0[2]} x_{1[0]} - x_{0[0]} x_{1[2]}
///      \\  x_{0[0]} x_{1[1]} - x_{0[1]} x_{1[0]}
///      \end{bmatrix}
/// \f]
///
template <typename Arith0, typename Arith1>
constexpr auto cross(
    const Array<Arith0, 3>& arr0, const Array<Arith1, 3>& arr1) noexcept {
    return Array<decltype(Arith0() * Arith1()), 3>{
        arr0[1] * arr1[2] - arr0[2] * arr1[1],
        arr0[2] * arr1[0] - arr0[0] * arr1[2],
        arr0[0] * arr1[1] - arr0[1] * arr1[0]};
}

/// 3-dimensional cross product.
///
/// \par Expression
/// \f[
///      \sum_{i,j,k} \varepsilon_{[i,j,k]} x_{0[i]} x_{1[j]} x_{2[k]}
///      \implies
///      \mathbf{x}_0 \cdot
///      (\mathbf{x}_1 \times \mathbf{x}_2)
/// \f]
///
template <typename Arith0, typename Arith1, typename Arith2>
constexpr auto cross(
    const Array<Arith0, 3>& arr0,
    const Array<Arith1, 3>& arr1,
    const Array<Arith2, 3>& arr2) noexcept {
    return dot(arr0, cross(arr1, arr2));
}

/// Euclidean length.
///
/// \par Expression
/// \f[
///      \lVert\mathbf{x}\rVert = \sqrt{\sum_k |x_{[k]}|^2}
/// \f]
///
/// \note
/// This pre-computes moduli to reduce everything to positive floating point
/// values. If squaring one of the reduced floating point values would cause
/// overflow or underflow, then this factors the maximum modulus out
/// from under the radical.
///
template <concepts::arithmetic_or_complex Arith, size_t N>
inline auto length(const Array<Arith, N>& arr) noexcept {
    if constexpr (N == 1) {
        return std::abs(arr[0]);
    }
    else if constexpr (N == 2) {
        return std::hypot(std::abs(arr[0]), std::abs(arr[1]));
    }
    else {
        using Float = to_floating_point_t<Arith>;
        Array<Float, N> tmp = pre::abs(arr);
        Float biggest = tmp[tmp.argmax()];
        if (biggest == 0)
            return Float(0);
        if (biggest * biggest >= Maximum<Float> / N ||
            biggest <= numeric_limits<Float>::min_squarable()) {
            if (biggest >= numeric_limits<Float>::min_invertible())
                tmp *= 1 / biggest;
            else
                tmp /= biggest; // Inverse overflows.
            return std::sqrt(dot(tmp, tmp)) * biggest;
        }
        else
            return std::sqrt(dot(tmp, tmp));
    }
}

/// Euclidean length squared.
///
/// \par Expression
/// \f[
///      \lVert\mathbf{x}\rVert^2 = \sum_k |x_{[k]}|^2
/// \f]
///
/// \note
/// This is implemented as `norm(arr).sum()`.
/// This is faster than but not as robust as `length()`, which accounts
/// for edge cases where entries overflow or underflow if squared. Since
/// the result here is length-squared anyway, additional precautions are
/// mostly unnecessary.
///
template <concepts::arithmetic_or_complex Arith, size_t N>
inline auto length2(const Array<Arith, N>& arr) noexcept {
    return pre::norm(arr).sum();
}

/// Euclidean distance.
template <concepts::arithmetic Arith0, concepts::arithmetic Arith1, size_t N>
inline auto distance(
    const Array<Arith0, N>& arr0, const Array<Arith1, N>& arr1) noexcept {
    return length(arr0 - arr1);
}

/// Euclidean distance squared.
template <concepts::arithmetic Arith0, concepts::arithmetic Arith1, size_t N>
inline auto distance2(
    const Array<Arith0, N>& arr0, const Array<Arith1, N>& arr1) noexcept {
    return length2(arr0 - arr1);
}

/// Normalize by Euclidean length.
template <concepts::arithmetic_or_complex Arith, size_t N>
inline auto normalize(const Array<Arith, N>& arr) noexcept {
    using Float = to_floating_point_t<Arith>;
    using Entry = decltype(arr[0] / Float(1));
    Array<Entry, N> res;
    if (Float len = length(arr); len > 0)
        res = arr / len;
    return res;
}

/// Normalize by Euclidean length, component version.
template <
    concepts::arithmetic_or_complex Arith,
    concepts::arithmetic_or_complex... Ariths>
inline auto normalize(Arith x0, Ariths... xs) noexcept {
    return normalize(Array{x0, xs...});
}

/// Fast (unsafe) normalize by Euclidean length.
template <concepts::arithmetic_or_complex Arith, size_t N>
inline auto fast_normalize(const Array<Arith, N>& arr) noexcept {
    return arr * (to_floating_point_t<Arith>(1) / pre::sqrt(length2(arr)));
}

/// Fast (unsafe) normalize by Euclidean length, component version.
template <
    concepts::arithmetic_or_complex Arith,
    concepts::arithmetic_or_complex... Ariths>
inline auto fast_normalize(Arith x0, Ariths... xs) noexcept {
    return fast_normalize(Array{x0, xs...});
}

/// Angle between vectors in 2-dimensions.
template <concepts::arithmetic Arith0, concepts::arithmetic Arith1>
inline auto angle_between(
    const Array<Arith0, 2>& arr0, const Array<Arith1, 2>& arr1) noexcept {
    return pre::atan2(cross(arr0, arr1), dot(arr0, arr1));
}

/// Matrix conjugate-transpose.
template <concepts::arithmetic_or_complex Arith, size_t M, size_t N>
constexpr Array<Arith, N, M> adjoint(const Array<Arith, M, N>& arr) noexcept {
    return pre::conj(transpose(arr));
}

/// Matrix determinant.
template <concepts::arithmetic_or_complex Arith, size_t N>
inline auto det(const Array<Arith, N, N>& arr) noexcept {
    if constexpr (N == 1) {
        return arr(0, 0);
    }
    else if constexpr (N == 2) {
        return cross(arr[0], arr[1]);
    }
    else if constexpr (N == 3) {
        return cross(arr[0], arr[1], arr[2]);
    }
    else {
        using Float = to_floating_point_t<Arith>;
        using Field = decltype(Arith() * Float());
        try {
            Array<Field, N, N> a = arr;
            Array<int, N> p;
            auto s = Linalg<Field>::lu(a, p); // Sign
            return a->diag().prod() * Field(s);
        }
        catch (const std::runtime_error&) {
            return Field(0);
        }
    }
}

/// Matrix inverse.
template <concepts::arithmetic_or_complex Arith, size_t N>
inline auto inverse(const Array<Arith, N, N>& arr) noexcept {
    if constexpr (std::integral<Arith>) {
        return inverse(Array<double, N, N>(arr));
    }
    else {
        using Float = to_floating_point_t<Arith>;
        using Field = Arith;
        if constexpr (N == 1) {
            return Array<Field, 1, 1>{Field(1) / arr[0][0]};
        }
        else if constexpr (N == 2) {
            Array<Field, 2, 2> cof = {
                +arr(1, 1), -arr(1, 0), //
                -arr(0, 1), +arr(0, 0)};
            return transpose(cof) * (Float(1) / dot(cof[0], arr[0]));
        }
        else if constexpr (N == 3) {
            Array<Field, 3, 3> cof;
            cof[0] = cross(arr[1], arr[2]);
            cof[1] = cross(arr[2], arr[0]);
            cof[2] = cross(arr[0], arr[1]);
            return transpose(cof) * (Float(1) / dot(cof[0], arr[0]));
        }
        else {
            try {
                Array<Field, N, N> a = arr;
                Array<Field, N, N> b = Array<Field, N, N>::identity();
                Array<Field, N, N> x;
                Array<int, N> p;
                Linalg<Field>::lu(a, p);
                Linalg<Field>::lu_solve(a, b, x, p);
                return x;
            }
            catch (const std::runtime_error&) {
                return Array<Field, N, N>(NaN<Float>);
            }
        }
    }
}

template <concepts::arithmetic_or_complex Arith, size_t M, size_t N>
inline auto svd(const Array<Arith, M, N>& arr) {
    using Float = to_floating_point_t<Arith>;
    using Field = decltype(Arith() * Float());
    Array<Field, M, N> x = arr;
    struct Result {
        Array<Float, pre::min(M, N)> s;
        Array<Field, M, M> u;
        Array<Field, N, N> v;
    } result;
    Linalg<Field> linalg;
    linalg.svd(x, result.u, result.v);
    result.s = Array<Field, pre::min(M, N)>(*x->diag());
    return result;
}

template <std::integral Int, size_t N>
constexpr Array<Int, N, N> permutation_matrix(
    const Array<Int, N>& p) noexcept {
    Array<Int, N, N> res;
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++)
            res(i, j) = p[i] == Int(j) ? 1 : 0;
    return res;
}

template <size_t N, std::floating_point Float>
constexpr Array<Float, N> linspace(Float a, Float b) noexcept {
    Array<Float, N> arr;
    for (size_t k = 0; k < N; k++)
        arr[k] = lerp(k / Float(N - 1), a, b);
    return arr;
    static_assert(N > 1);
}

/// Initializers for 2-dimensional floating point arrays.
///
/// \note
/// To access these initializers, use `pre::Vec2<Float>::`.
///
template <std::floating_point Float>
struct Array_initializers<Array<Float, 2>> {
    /// Uniform disk probability density function.
    static constexpr Float uniform_disk_pdf() noexcept {
        return M_1_pi<Float>;
    }

    /// Uniform disk probability density function sampling routine.
    static Vec2<Float> uniform_disk_pdf_sample(Vec2<Float> u) noexcept {
        u = 2 * u - 1;
        if ((u == 0).all()) {
            return u;
        }
        else {
            Float r;
            Float theta;
            if (pre::abs(u[0]) > pre::abs(u[1])) {
                r = u[0];
                theta = M_pi_4<Float> * (u[1] / u[0]);
            }
            else {
                r = u[1];
                theta = M_pi_4<Float> * (u[0] / u[1]);
                theta = M_pi_2<Float> - theta;
            }
            return {r * pre::cos(theta), r * pre::sin(theta)};
        }
    }
};

// TODO Move this
/// Initializers for 3-dimensional floating point arrays.
///
/// \note
/// To access these initializers, use `pre::Vec3<Float>::`.
///
template <std::floating_point Float>
struct Array_initializers<Array<Float, 3>> {

    /// Form spherical direction from z-component and azimuth angle.
    static Vec3<Float> sphericalz(Float z, Float phi) noexcept {
        Float cos_theta = clamp_abs(z, +1);
        Float sin_theta = pre::sqrt(1 - z * z);
        auto [sin_phi, cos_phi] = sincos(phi);
        return {sin_theta * cos_phi, sin_theta * sin_phi, cos_theta};
    }

    /// Uniform hemisphere probability density function.
    static constexpr Float uniform_hemisphere_pdf() noexcept {
        return M_1_pi<Float> / 2;
    }

    /// Uniform hemisphere probability density function sampling routine.
    static Vec3<Float> uniform_hemisphere_pdf_sample(Vec2<Float> u) noexcept {
        return sphericalz(u[0], 2 * M_pi<Float> * u[1]);
    }

    /// Uniform sphere probability density function.
    static constexpr Float uniform_sphere_pdf() noexcept {
        return M_1_pi<Float> / 4;
    }

    /// Uniform sphere probability density function sampling routine.
    static Vec3<Float> uniform_sphere_pdf_sample(Vec2<Float> u) noexcept {
        return sphericalz(2 * u[0] - 1, 2 * M_pi<Float> * u[1]);
    }

    /// Cosine hemisphere probability density function.
    static constexpr Float cosine_hemisphere_pdf(Float z) noexcept {
        return M_1_pi<Float> * pre::max(z, 0);
    }

    /// Cosine hemisphere probability density function sampling routine.
    static Vec3<Float> cosine_hemisphere_pdf_sample(Vec2<Float> u) noexcept {
        Vec2<Float> p = Vec2<Float>::uniform_disk_pdf_sample(u);
        Vec3<Float> w = {p[0], p[1], pre::sqrt(1 - pre::min(length2(p), 1))};
        return w;
    }

    /// Uniform cone probability density function.
    static constexpr Float uniform_cone_pdf(Float zmax) noexcept {
        return (M_1_pi<Float> / 2) / (1 - zmax);
    }

    /// Uniform cone probability density function sampling routine.
    static Vec3<Float> uniform_cone_pdf_sample(
        Float zmax, Vec2<Float> u) noexcept {
        return sphericalz((1 - u[0]) + u[0] * zmax, 2 * M_pi<Float> * u[1]);
    }

    /// Henyey-Greenstein phase probability density function.
    static Float hg_phase_pdf(Float g, Float z) noexcept {
        if (pre::abs(g) < Float(0.00001)) {
            return uniform_sphere_pdf();
        }
        else {
            g = clamp_abs(g, 0.99999);
            Float a = 1 - g * g;
            Float b = 1 + g * g - 2 * g * z;
            return (M_1_pi<Float> / 4) * (a / (b * pre::sqrt(b)));
        }
    }

    /// Henyey-Greenstein phase probability density function sampling routine.
    static Vec3<Float> hg_phase_pdf_sample(Float g, Vec2<Float> u) noexcept {
        if (pre::abs(g) < Float(0.00001)) {
            return uniform_sphere_pdf_sample(u);
        }
        else {
            g = clamp_abs(g, 0.99999);
            Float h = (1 - g * g) / (1 - g + 2 * g * u[0]);
            Float z = (1 + g * g - h * h) / (2 * g);
            return sphericalz(z, 2 * M_pi<Float> * u[1]);
        }
    }
};

/// Initializers for 2x2-dimensional floating point arrays.
///
/// \note
/// To access these initializers, use `pre::Mat2<Float>::`.
///
template <std::floating_point Float>
struct Array_initializers<Array<Float, 2, 2>> {
    /// Rotate counter-clockwise.
    static Mat2<Float> rotate(Float theta) noexcept {
        auto [sin_theta, cos_theta] = sincos(theta);
        return {
            {+cos_theta, -sin_theta}, //
            {+sin_theta, +cos_theta}};
    }
};

/// Initializers for 3x3-dimensional floating point arrays.
///
/// \note
/// To access these initializers, use `pre::Mat3<Float>::`.
///
template <std::floating_point Float>
struct Array_initializers<Array<Float, 3, 3>> {
    /// Build real 3-dimensional orthonormal basis.
    ///
    /// \par Expression
    /// - \f$ \alpha_0 \gets -1 / (\hat{z}_{[2]} + 1) \f$
    /// - \f$ \alpha_1 \gets \alpha_0 \hat{z}_{[0]} \hat{z}_{[1]} \f$
    /// - \f$ \alpha_2 \gets \alpha_0 \hat{z}_{[0]}^2 + 1 \f$
    /// - \f$ \alpha_3 \gets \alpha_0 \hat{z}_{[1]}^2 + 1 \f$
    /// - \f$ \hat{\mathbf{x}} \gets
    ///   [\alpha_2\; \alpha_1\; -\hat{z}_{[0]}]^\top \f$
    /// - \f$ \hat{\mathbf{y}} \gets
    ///   [\alpha_1\; \alpha_3\; -\hat{z}_{[1]}]^\top \f$
    ///
    /// \note
    /// As the notation suggests, the implementation assumes the input
    /// vector `hatz` is unit-length.
    ///
    static Mat3<Float> build_onb(Vec3<Float> hatz) noexcept {
        Vec3<Float> hatx = {};
        Vec3<Float> haty = {};
        if (hatz[2] < Float(-0.9999999)) {
            hatz[2] = -1;
            hatx[1] = -1;
            haty[0] = -1;
        }
        else {
            Float alpha0 = -1 / (hatz[2] + 1);
            Float alpha1 = alpha0 * hatz[0] * hatz[1];
            Float alpha2 = alpha0 * hatz[0] * hatz[0] + 1;
            Float alpha3 = alpha0 * hatz[1] * hatz[1] + 1;
            hatx = {alpha2, alpha1, -hatz[0]};
            haty = {alpha1, alpha3, -hatz[1]};
            if (hatz[2] < Float(-0.999)) {
                // Accuracy begins to suffer, so manually orthonormalize.
                hatx = normalize(hatx - dot(hatx, hatz) * hatz);
                haty = normalize(
                    haty - dot(haty, hatx) * hatx - dot(haty, hatz) * hatz);
            }
        }
        return {
            {hatx[0], haty[0], hatz[0]}, //
            {hatx[1], haty[1], hatz[1]}, //
            {hatx[2], haty[2], hatz[2]}};
    }

    /// Rotate counter-clockwise around arbitrary axis.
    ///
    /// \param[in] theta  Angle in radians.
    /// \param[in] hatv   Normalized rotation axis.
    ///
    /// \par Expression
    /// \f[
    ///      \begin{bmatrix}
    ///          v_x^2   (1 - \cos{\theta}) + \cos{\theta}
    ///      &   v_x v_y (1 - \cos{\theta}) - v_z \sin{\theta}
    ///      &   v_x v_z (1 - \cos{\theta}) + v_y \sin{\theta}
    ///      \\  v_x v_y (1 - \cos{\theta}) + v_z \cos{\theta}
    ///      &   v_y^2   (1 - \cos{\theta}) + \cos{\theta}
    ///      &   v_y v_z (1 - \cos{\theta}) - v_x \sin{\theta}
    ///      \\  v_x v_z (1 - \cos{\theta}) - v_y \sin{\theta}
    ///      &   v_y v_z (1 - \cos{\theta}) + v_x \sin{\theta}
    ///      &   v_z^2   (1 - \cos{\theta}) + \cos{\theta}
    ///      \end{bmatrix}
    /// \f]
    ///
    static Mat3<Float> rotate(Float theta, Vec3<Float> hatv) noexcept {
        auto [sin_theta, cos_theta] = sincos(theta);
        Float vx = hatv[0];
        Float vy = hatv[1];
        Float vz = hatv[2];
        Float vxvx = vx * vx, vxvy = vx * vy, vxvz = vx * vz;
        Float vyvy = vy * vy, vyvz = vy * vz;
        Float vzvz = vz * vz;
        return {
            {vxvx * (1 - cos_theta) + cos_theta,
             vxvy * (1 - cos_theta) - vz * sin_theta,
             vxvz * (1 - cos_theta) + vy * sin_theta}, //
            {vxvy * (1 - cos_theta) + vz * sin_theta,
             vyvy * (1 - cos_theta) + cos_theta,
             vyvz * (1 - cos_theta) - vx * sin_theta}, //
            {vxvz * (1 - cos_theta) - vy * sin_theta,
             vyvz * (1 - cos_theta) + vx * sin_theta,
             vzvz * (1 - cos_theta) + cos_theta}};
    }

    /// Rotate counter-clockwise around X-axis.
    static Mat3<Float> rotatex(Float theta) noexcept {
        return rotate(theta, {1, 0, 0});
    }

    /// Rotate counter-clockwise around Y-axis.
    static Mat3<Float> rotatey(Float theta) noexcept {
        return rotate(theta, {0, 1, 0});
    }

    /// Rotate counter-clockwise around Z-axis.
    static Mat3<Float> rotatez(Float theta) noexcept {
        return rotate(theta, {0, 0, 1});
    }

    /// Uniform scale.
    static Mat3<Float> scale(Float alpha) noexcept {
        return {{alpha, 0, 0}, {0, alpha, 0}, {0, 0, alpha}};
    }

    /// Non-uniform scale along arbitrary axis.
    static Mat3<Float> scale(Float alpha, Vec3<Float> hatv) noexcept {
        Mat3<Float> res = (alpha - 1) * outer(hatv, hatv);
        res->diag() += 1;
        return res;
    }

    /// Non-uniform axis-aligned scale.
    static Mat3<Float> scale(Vec3<Float> alpha) noexcept {
        return {
            {alpha[0], 0, 0}, //
            {0, alpha[1], 0}, //
            {0, 0, alpha[2]}};
    }
};

/// Initializers for 4x4-dimensional floating point arrays.
///
/// \note
/// To access these initializers, use `pre::Mat4<Float>::`.
///
template <std::floating_point Float>
struct Array_initializers<Array<Float, 4, 4>> {
    /// Translate.
    static Mat4<Float> translate(Vec3<Float> v) noexcept {
        return {
            {1, 0, 0, v[0]}, //
            {0, 1, 0, v[1]}, //
            {0, 0, 1, v[2]}, //
            {0, 0, 0, 1}};
    }

    /// Rotate counter-clockwise around arbitrary axis.
    static Mat4<Float> rotate(Float theta, Vec3<Float> hatv) noexcept {
        Mat4<Float> res;
        res(Slice(0, 3), Slice(0, 3)) = **Mat3<Float>::rotate(theta, hatv);
        res(3, 3) = 1;
        return res;
    }

    /// Rotate counter-clockwise around X-axis.
    static Mat4<Float> rotatex(Float theta) noexcept {
        return rotate(theta, {1, 0, 0});
    }

    /// Rotate counter-clockwise around Y-axis.
    static Mat4<Float> rotatey(Float theta) noexcept {
        return rotate(theta, {0, 1, 0});
    }

    /// Rotate counter-clockwise around Z-axis.
    static Mat4<Float> rotatez(Float theta) noexcept {
        return rotate(theta, {0, 0, 1});
    }

    /// Look-at transform.
    ///
    /// \param[in] pfrom  Look-from point.
    /// \param[in] pto    Look-to point.
    /// \param[in] vup    Up vector.
    ///
    /// \par Expression
    /// \f[
    ///      \begin{bmatrix}
    ///          \hat{\mathbf{x}}^\top
    ///      &  -\hat{\mathbf{x}}^\top \mathbf{p}_{\text{from}}
    ///      \\  \hat{\mathbf{y}}^\top
    ///      &  -\hat{\mathbf{y}}^\top \mathbf{p}_{\text{from}}
    ///      \\  \hat{\mathbf{z}}^\top
    ///      &  -\hat{\mathbf{z}}^\top \mathbf{p}_{\text{from}}
    ///      \\  \mathbf{0}^\top & 1
    ///      \end{bmatrix}
    /// \f]
    /// where
    /// - \f$\mathbf{z}\gets\mathbf{p}_{\text{from}}-\mathbf{p}_{\text{to}}\f$
    /// - \f$\mathbf{x}\gets\mathbf{v}_{\text{up}}\times\mathbf{z}\f$
    /// - \f$\hat{\mathbf{z}}\gets\operatorname{normalize}(\mathbf{z})\f$
    /// - \f$\hat{\mathbf{x}}\gets\operatorname{normalize}(\mathbf{x})\f$
    /// - \f$\hat{\mathbf{y}}\gets\hat{\mathbf{z}}\times\hat{\mathbf{x}}\f$
    ///
    /// \note
    /// This implementation is consistent with OpenGL conventions.
    ///
    static Mat4<Float> look_at(
        const Vec3<Float>& pfrom,
        const Vec3<Float>& pto,
        const Vec3<Float>& vup) noexcept {
        Vec3<Float> z = pfrom - pto;
        Vec3<Float> x = cross(vup, z);
        Vec3<Float> hatz = normalize(z);
        Vec3<Float> hatx = normalize(x);
        Vec3<Float> haty = cross(hatz, hatx);
        return {
            {hatx[0], hatx[1], hatx[2], -dot(hatx, pfrom)},
            {haty[0], haty[1], haty[2], -dot(haty, pfrom)},
            {hatz[0], hatz[1], hatz[2], -dot(hatz, pfrom)},
            {0, 0, 0, 1}};
    }

    /// Orthographic projection.
    ///
    /// \param[in] x0  Frustum X-coordinate 0.
    /// \param[in] x1  Frustum X-coordinate 1.
    /// \param[in] y0  Frustum Y-coordinate 0.
    /// \param[in] y1  Frustum Y-coordinate 1.
    /// \param[in] z0  Frustum Z-coordinate 0.
    /// \param[in] z1  Frustum Z-coordinate 1.
    ///
    /// \par Expression
    /// \f[
    ///      \begin{bmatrix}
    ///          2 / (x1 - x0) & 0 & 0 & -(x1 + x0) / (x1 - x0)
    ///      \\  0 & 2 / (y1 - y0) & 0 & -(y1 + y0) / (y1 - y0)
    ///      \\  0 & 0 &-2 / (z1 - z0) & -(z1 + z0) / (z1 - z0)
    ///      \\  0 & 0 & 0 & 1
    ///      \end{bmatrix}
    /// \f]
    ///
    /// \note
    /// This implementation is consistent with OpenGL conventions.
    ///
    static Mat4<Float> ortho(
        Float x0, Float x1, Float y0, Float y1, Float z0, Float z1) noexcept {
        Mat4<Float> res;
        res->diag() = {2 / (x1 - x0), 2 / (y1 - y0), -2 / (z1 - z0), 1};
        res(Slice(0, 3), 3) = {
            -(x1 + x0) / (x1 - x0), //
            -(y1 + y0) / (y1 - y0), //
            -(z1 + z0) / (z1 - z0)};
        return res;
    }

    /// Perspective projection.
    ///
    /// \param[in] x0  Frustum X-coordinate 0 at Z-coordinate 0.
    /// \param[in] x1  Frustum X-coordinate 1 at Z-coordinate 0.
    /// \param[in] y0  Frustum Y-coordinate 0 at Z-coordinate 0.
    /// \param[in] y1  Frustum Y-coordinate 1 at Z-coordinate 0.
    /// \param[in] z0  Frustum Z-coordinate 0.
    /// \param[in] z1  Frustum Z-coordinate 1.
    ///
    /// \par Expression
    /// \f[
    ///      \begin{bmatrix}
    ///          2 z_0 / (x_1 - x_0) & 0 & (x_1 + x_0) / (x_1 - x_0) & 0
    ///      \\  0 & 2 z_0 / (y_1 - y_0) & (y_1 + y_0) / (y_1 - y_0) & 0
    ///      \\  0 & 0 & -(z_1 + z_0) / (z_1 - z_0) & -2 z_1 z_0 / (z_1 - z_0)
    ///      \\  0 & 0 & -1 & 0
    ///      \end{bmatrix}
    /// \f]
    ///
    /// \note
    /// This implementation is consistent with OpenGL conventions.
    ///
    static Mat4<Float> persp(
        Float x0, Float x1, Float y0, Float y1, Float z0, Float z1) noexcept {
        Mat4<Float> res;
        res->diag() = {
            2 * z0 / (x1 - x0), //
            2 * z0 / (y1 - y0), //
            -(z1 + z0) / (z1 - z0), 0};
        res(Slice(0, 2), 2) = {(x1 + x0) / (x1 - x0), (y1 + y0) / (y1 - y0)};
        res(2, 3) = -2 * z1 * z0 / (z1 - z0);
        res(3, 2) = -1;
        return res;
    }

    /// Perspective projection.
    ///
    /// \param[in] f   Frustum Y-FOV in radians.
    /// \param[in] r   Frustum X-to-Y ratio.
    /// \param[in] z0  Frustum Z-coordinate 0.
    /// \param[in] z1  Frustum Z-coordinate 1.
    ///
    /// \par Expression
    /// \f[
    ///      \begin{bmatrix}
    ///          \cot(f / 2) / r & 0 & 0 & 0
    ///      \\  0 & \cot(f / 2) & 0 & 0
    ///      \\  0 & 0 & -(z_1 + z_0) / (z_1 - z_0) & -2 z_1 z_0 / (z_1 - z_0)
    ///      \\  0 & 0 & -1 & 0
    ///      \end{bmatrix}
    /// \f]
    ///
    static Mat4<Float> persp(Float f, Float r, Float z0, Float z1) noexcept {
        Mat4<Float> res;
        Float cot_f_2 = cot(f / 2);
        res->diag() = {cot_f_2 / r, cot_f_2, -(z1 + z0) / (z1 - z0), 0};
        res(2, 3) = -2 * z1 * z0 / (z1 - z0);
        res(3, 2) = -1;
        return res;
    }
};

} // namespace pre
