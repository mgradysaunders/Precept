/*-*- C++ -*-*/
#pragma once

namespace pre {

/// Numeric constants, empty by default.
template <typename>
struct numeric_constants {};

/// Numeric constants, `std::floating_point` specialization.
template <std::floating_point Float>
struct numeric_constants<Float> {
    /// \f$ e \f$.
    static consteval Float M_e() noexcept {
        return Float(2.7182818284590452353602874713526625L);
    }

    /// \f$ \log_2(e) \f$.
    static consteval Float M_log2e() noexcept {
        return Float(1.4426950408889634073599246810018921L);
    }

    /// \f$ \log_{10}(e) \f$.
    static consteval Float M_log10e() noexcept {
        return Float(0.4342944819032518276511289189166051L);
    }

    /// \f$ \log_e(2) \f$.
    static consteval Float M_ln2() noexcept {
        return Float(0.6931471805599453094172321214581766L);
    }

    /// \f$ \log_e(10) \f$.
    static consteval Float M_ln10() noexcept {
        return Float(2.3025850929940456840179914546843642L);
    }

    /// \f$ \pi \f$.
    static consteval Float M_pi() noexcept {
        return Float(3.1415926535897932384626433832795029L);
    }

    /// \f$ \pi/2 \f$.
    static consteval Float M_pi_2() noexcept {
        return Float(1.5707963267948966192313216916397514L);
    }

    /// \f$ \pi/4 \f$.
    static consteval Float M_pi_4() noexcept {
        return Float(0.7853981633974483096156608458198757L);
    }

    /// \f$ 1/\pi \f$.
    static consteval Float M_1_pi() noexcept {
        return Float(0.3183098861837906715377675267450287L);
    }

    /// \f$ 2/\pi \f$.
    static consteval Float M_2_pi() noexcept {
        return Float(0.6366197723675813430755350534900574L);
    }

    /// \f$ 1/\sqrt{\pi} \f$.
    static consteval Float M_1_sqrtpi() noexcept {
        return Float(1.1283791670955125738961589031215452L) * Float(0.5);
    }

    /// \f$ 2/\sqrt{\pi} \f$.
    static consteval Float M_2_sqrtpi() noexcept {
        return Float(1.1283791670955125738961589031215452L);
    }

    /// \f$ \sqrt{2} \f$.
    static consteval Float M_sqrt2() noexcept {
        return Float(1.4142135623730950488016887242096981L);
    }

    /// \f$ \sqrt{1/2} \f$.
    static consteval Float M_sqrt1_2() noexcept {
        return Float(0.7071067811865475244008443621048490L);
    }

    /// \f$ \gamma \f$ (Euler's constant).
    ///
    /// \f[
    ///     \gamma =
    ///         \lim_{n\to\infty}
    ///         \left(-\log(n) + \sum_{k=1}^{n}\frac{1}{k}\right)
    /// \f]
    ///
    static consteval Float M_gamma() noexcept {
        return Float(0.5772156649015328606065120900824024L);
    }

    /// \f$ h \f$ (Planck's constant).
    ///
    /// \f[
    ///     h = 6.62607015\times10^{-34}\,\mathrm{J}\cdot\mathrm{s}
    /// \f]
    //
    static consteval Float M_h() noexcept {
        return Float(6.62607015e-34L);
    }

    ///
    /// \f$ c \f$ (light speed).
    ///
    /// \f[
    ///     c = 299792458\,\mathrm{m}/\mathrm{s}
    /// \f]
    ///
    static consteval Float M_c() noexcept {
        return Float(299792458);
    }
};

template <std::floating_point Float>
constexpr Float M_pi = numeric_constants<Float>::M_pi();

template <std::floating_point Float>
constexpr Float M_pi_2 = numeric_constants<Float>::M_pi_2();

template <std::floating_point Float>
constexpr Float M_pi_4 = numeric_constants<Float>::M_pi_4();

template <std::floating_point Float>
constexpr Float M_1_pi = numeric_constants<Float>::M_1_pi();

template <std::floating_point Float>
constexpr Float M_2_pi = numeric_constants<Float>::M_2_pi();

template <std::floating_point Float>
constexpr Float M_1_sqrtpi = numeric_constants<Float>::M_1_sqrtpi();

template <std::floating_point Float>
constexpr Float M_2_sqrtpi = numeric_constants<Float>::M_2_sqrtpi();

template <std::floating_point Float>
constexpr Float M_sqrt2 = numeric_constants<Float>::M_sqrt2();

template <std::floating_point Float>
constexpr Float M_sqrt1_2 = numeric_constants<Float>::M_sqrt1_2();

/// Numeric constants, `std::complex` specialization.
template <std::floating_point Float>
struct numeric_constants<std::complex<Float>> : numeric_constants<Float> {};

} // namespace pre

constexpr auto operator""_degreesf(long double x) noexcept {
    return float(x * pre::numeric_constants<long double>::M_pi() / 180.0L);
}

constexpr auto operator""_degrees(long double x) noexcept {
    return double(x * pre::numeric_constants<long double>::M_pi() / 180.0L);
}
