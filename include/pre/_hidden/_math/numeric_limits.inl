/*-*- C++ -*-*/
#pragma once

namespace pre {

/// Wrap `std::numeric_limits`.
template <typename T>
struct numeric_limits : std::numeric_limits<T> {};

/// Wrap `std::numeric_limits`, with extra floating point constants.
template <std::floating_point T>
struct numeric_limits<T> : std::numeric_limits<T> {
    /// Machine epsilon.
    ///
    /// This is the greatest positive floating point value \f$ x \f$
    /// that satisfies \f$ 1 \oplus x = 1 \f$.
    ///
    static consteval T machine_epsilon() noexcept {
        // Correct for IEEE floating point types, possibly others?
        return numeric_limits::epsilon() / 2;
    }

    /// Echelon.
    ///
    /// \par Expression
    /// \f[
    ///     \tilde{\epsilon}(n) =
    ///         \frac{\epsilon_{\text{m}} n}{1 -
    ///               \epsilon_{\text{m}} n}
    /// \f]
    /// where \f$ \epsilon_{\text{m}} \f$ is machine epsilon.
    ///
    static consteval T echelon(unsigned int n) noexcept {
        return machine_epsilon() * n / (1 - machine_epsilon() * n);
    }

    /// Minimum squarable value.
    ///
    /// This is the least positive floating point value \f$ x \f$
    /// that satisfies \f$ x \otimes x > 0 \f$.
    ///
    static consteval T min_squarable() noexcept {
        if constexpr (std::is_same_v<T, float>) {
            return 0x1.000002p-75f; // 2^(-75)
        }
        else if constexpr (std::is_same_v<T, double>) {
            return 0x1.6a09e667f3bcdp-538; // 2^(-537.5)
        }
        else if constexpr (std::is_same_v<T, long double>) {
            return 0x8.000000000000001p-8226L; // 2^(-8223)
        }
        else {
            // Error.
            throw std::logic_error("unsupported floating point type");
        }
    }

    /// Minimum invertible value.
    ///
    /// This is the least positive floating point value \f$ x \f$
    /// that satisfies \f$ 1 \oslash x < \infty \f$.
    ///
    static consteval T min_invertible() noexcept {
        // Correct for IEEE floating point types, possibly others?
        return numeric_limits::min() / 4 + numeric_limits::denorm_min();
    }
};

/// Wrap `std::numeric_limits`, `std::complex` specialization.
template <std::floating_point T>
struct numeric_limits<std::complex<T>> : numeric_limits<T> {};

template <typename T>
constexpr auto Minimum = numeric_limits<T>::min();

template <typename T>
constexpr auto Maximum = numeric_limits<T>::max();

template <typename T>
constexpr auto Inf = numeric_limits<T>::infinity();

template <typename T>
constexpr auto NaN = numeric_limits<T>::quiet_NaN();

template <typename T>
constexpr auto Eps = numeric_limits<T>::epsilon();

template <typename T>
constexpr auto MachineEps = numeric_limits<T>::machine_epsilon();

template <typename T, unsigned n>
constexpr auto MachineEch = numeric_limits<T>::echelon(n);

} // namespace pre
