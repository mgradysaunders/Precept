/*-*- C++ -*-*/
#pragma once

namespace pre {

/// \name Common
/** \{ */

/// Analogous to `std::min()`, except only for numbers.
template <
    concepts::arithmetic_or_enum Arith0,
    concepts::arithmetic_or_enum Arith1>
constexpr auto min(Arith0 x, Arith1 y) noexcept {
    return x < y ? x : y;
}

/// Analogous to `std::max()`, except only for numbers.
template <
    concepts::arithmetic_or_enum Arith0,
    concepts::arithmetic_or_enum Arith1>
constexpr auto max(Arith0 x, Arith1 y) noexcept {
    return x < y ? y : x;
}

/// Equivalent to `x = min(x, y)`.
template <
    concepts::arithmetic_or_enum Arith0,
    concepts::arithmetic_or_enum Arith1>
constexpr void minimize(Arith0& x, Arith1 y) noexcept {
    x = pre::min(x, y);
}

/// Equivalent to `x = max(x, y)`.
template <
    concepts::arithmetic_or_enum Arith0,
    concepts::arithmetic_or_enum Arith1>
constexpr void maximize(Arith0& x, Arith1 y) noexcept {
    x = pre::max(x, y);
}

/// Clamp in range.
///
/// \param[in] x  Value.
/// \param[in] a  Range minimum.
/// \param[in] b  Range maximum.
///
template <concepts::arithmetic_or_enum Arith>
constexpr Arith clamp(
    Arith x,
    std::type_identity_t<Arith> a,
    std::type_identity_t<Arith> b) noexcept {
    return pre::min(pre::max(x, a), b);
}

/// Clamp absolute value in range.
///
/// \param[in] x  Value.
/// \param[in] b  Range maximum.
///
template <concepts::arithmetic_or_enum Arith>
constexpr Arith clamp_abs(Arith x, std::type_identity_t<Arith> b) noexcept {
    return clamp<Arith>(x, -b, +b);
}

/// Raise to integer power.
template <concepts::arithmetic Arith>
constexpr Arith nthpow(Arith x, int n) noexcept {
    if (n < 0) {
        return Arith(1) / nthpow(x, -n);
    }
    else {
        Arith y = 1;
        while (1) {
            if (n & 1)
                y *= x;
            n >>= 1;
            if (n == 0)
                break;
            x *= x;
        }
        return y;
    }
}

/** \} */

/// \name Integer helpers
/** \{ */

/// Wrap integer in range.
///
/// - If \f$ n > 0 \f$, wrap \f$ k \f$ to \f$ [0, n) \f$.
/// - If \f$ n < 0 \f$, wrap \f$ k \f$ to \f$ (n, 0] \f$.
///
/// \note
/// For \f$ k > 0 \f$, \f$ n > 0 \f$, this operation is
/// equivalent to \f$ k \% n \f$.
///
template <std::integral Int>
constexpr Int repeat(Int k, Int n) noexcept {
    if constexpr (std::unsigned_integral<Int>) {
        return k % n;
    }
    else {
        if (n < Int(0))
            return -repeat(-k, -n);
        if (k >= Int(0))
            return k % n;
        k = n + k % n;
        if (k == n)
            k = Int(0);
        return k;
    }
}

/// Wrap integer in range and mirror with each repeat.
///
/// - If \f$ n > 0 \f$, wrap \f$ k \f$ to \f$ [0, n) \f$.
/// - If \f$ n < 0 \f$, wrap \f$ k \f$ to \f$ (n, 0] \f$.
///
template <std::integral Int>
constexpr Int mirror(Int k, Int n) noexcept {
    if (n < Int(0)) {
        return -mirror(-k, -n);
    }
    else {
        Int r = k % n;
        Int q = k / n;
        if (r < Int(0)) {
            r += n;
            q++;
        }
        if (q & Int(1)) {
            r = n - r - Int(1);
        }
        return r;
    }
}

/// Is power of 2?
template <std::integral Int>
constexpr bool ispow2(Int n) noexcept {
    return n > 0 && (n & (n - 1)) == 0;
}

/// Round up to power of 2.
template <std::integral Int>
constexpr Int roundpow2(Int n) noexcept {
    if (n <= 0) {
        return 1;
    }
    else {
        n--;
        Int p = 1;
        while (n) {
            n >>= 1;
            p <<= 1;
        }
        return p;
    }
}

/// Find first bit equal to 1.
///
/// \note
/// This is equal to the binary logarithm if the input integer
/// is a positive power of 2.
///
/// \note
/// Uses [`__builtin_ctz`][1] when compiling with GCC.
/// [1]: https://gcc.gnu.org/onlinedocs/gcc/Other-Builtins.html
///
template <std::integral Int>
constexpr Int first1(Int n) noexcept {
    if (n == 0)
        return 0; // Error?
#if __GNUC__
    if constexpr (sizeof(Int) <= sizeof(int))
        return __builtin_ctz(n);
    else if constexpr (sizeof(Int) == sizeof(long))
        return __builtin_ctzl(n);
    else if constexpr (sizeof(Int) == sizeof(long long))
        return __builtin_ctzll(n);
#endif // #if __GNUC__

    Int j = 0;
    while (!(n & 1)) {
        n >>= 1;
        j++;
    }
    return j;
}

/// Cyclical bit rotate left.
template <std::integral Int>
constexpr Int rotl(Int val, unsigned rot) noexcept {
    return (val << rot) | (val >> ((-rot) & (sizeof(Int) * 8 - 1)));
}

/// Cyclical bit rotate right.
template <std::integral Int>
constexpr Int rotr(Int val, unsigned rot) noexcept {
    return (val >> rot) | (val << ((-rot) & (sizeof(Int) * 8 - 1)));
}

/// Bit swap.
template <std::integral Int>
constexpr Int bit_swap(Int val, int pos0, int pos1) noexcept {
    Int mask = ((val >> pos0) ^ (val >> pos1)) & 1;
    return ((mask << pos0) | (mask << pos1)) ^ val;
}

/// Bit reverse.
template <std::integral Int>
constexpr Int bit_reverse(Int val) noexcept {
    if constexpr (std::signed_integral<Int>) {
        return bit_reverse<std::make_unsigned_t<Int>>(val);
    }
    else if constexpr (sizeof(Int) == sizeof(std::uint8_t)) {
        constexpr std::uint8_t mask[3] = {0xaaU, 0xccU, 0xf0U};
        for (int k = 0; k < 3; k++)
            val =
                ((val & mask[k]) >> (1 << k)) | ((val & ~mask[k]) << (1 << k));
    }
    else if constexpr (sizeof(Int) == sizeof(std::uint16_t)) {
        constexpr std::uint16_t mask[4] = {0xaaaaU, 0xccccU, 0xf0f0U, 0xff00U};
#if __GNUC__
        val = __builtin_bswap16(val);
        for (int k = 0; k < 3; k++)
#else
        for (int k = 0; k < 4; k++)
#endif // #if __GNUC__
            val =
                ((val & mask[k]) >> (1 << k)) | ((val & ~mask[k]) << (1 << k));
    }
    else if constexpr (sizeof(Int) == sizeof(std::uint32_t)) {
        constexpr std::uint32_t mask[5] = {
            0xaaaaaaaaUL, 0xccccccccUL, 0xf0f0f0f0UL, 0xff00ff00UL,
            0xffff0000UL};
#if __GNUC__
        val = __builtin_bswap32(val);
        for (int k = 0; k < 3; k++)
#else
        for (int k = 0; k < 5; k++)
#endif // #if __GNUC__
            val =
                ((val & mask[k]) >> (1 << k)) | ((val & ~mask[k]) << (1 << k));
    }
    else if constexpr (sizeof(Int) == sizeof(std::uint64_t)) {
        constexpr std::uint64_t mask[6] = {
            0xaaaaaaaaaaaaaaaaULL, 0xccccccccccccccccULL,
            0xf0f0f0f0f0f0f0f0ULL, 0xff00ff00ff00ff00ULL,
            0xffff0000ffff0000ULL, 0xffffffff00000000ULL};
#if __GNUC__
        val = __builtin_bswap64(val);
        for (int k = 0; k < 3; k++)
#else
        for (int k = 0; k < 6; k++)
#endif // #if __GNUC__
            val =
                ((val & mask[k]) >> (1 << k)) | ((val & ~mask[k]) << (1 << k));
    }
    return val;
}

// Bit interleave with zero.
template <std::integral Int>
constexpr Int bit_interleave_zero(Int val) noexcept {
    if constexpr (std::signed_integral<Int>) {
        return bit_interleave_zero<std::make_unsigned_t<Int>>(val);
    }
    else if constexpr (sizeof(Int) == sizeof(std::uint8_t)) {
        constexpr std::uint8_t mask[3] = {0xaaU, 0xccU, 0xf0U};
        for (int k = 2; k >= 0; k--)
            val = (val ^ (val << (1 << k))) & ~mask[k];
    }
    else if constexpr (sizeof(Int) == sizeof(std::uint16_t)) {
        constexpr std::uint16_t mask[4] = {0xaaaaU, 0xccccU, 0xf0f0U, 0xff00U};
        for (int k = 3; k >= 0; k--)
            val = (val ^ (val << (1 << k))) & ~mask[k];
    }
    else if constexpr (sizeof(Int) == sizeof(std::uint32_t)) {
        constexpr std::uint32_t mask[5] = {
            0xaaaaaaaaUL, 0xccccccccUL, //
            0xf0f0f0f0UL, 0xff00ff00UL, //
            0xffff0000UL};
        for (int k = 4; k >= 0; k--)
            val = (val ^ (val << (1 << k))) & ~mask[k];
    }
    else if constexpr (sizeof(Int) == sizeof(std::uint64_t)) {
        constexpr std::uint64_t mask[6] = {
            0xaaaaaaaaaaaaaaaaULL, 0xccccccccccccccccULL,
            0xf0f0f0f0f0f0f0f0ULL, 0xff00ff00ff00ff00ULL,
            0xffff0000ffff0000ULL};
        for (int k = 5; k >= 0; k--)
            val = (val ^ (val << (1 << k))) & ~mask[k];
    }
    return val;
}

/// Bit interleave.
///
/// Form an integer \f$ z \f$ from the
/// bits of \f$ x \f$ and \f$ y \f$ such that
/// - bit \f$ k \f$ of \f$ x \f$ is bit \f$ 2k \f$ of \f$ z \f$,
/// - bit \f$ k \f$ of \f$ y \f$ is bit \f$ 2k + 1 \f$ of \f$ z \f$.
///
/// \note
/// If the input integer type has \f$ 2^n \f$ bits, the top
/// \f$ 2^{n - 1} \f$ bits of the input integers do not appear in
/// the output.
///
template <std::integral Int>
constexpr Int bit_interleave(Int val0, Int val1) noexcept {
    val0 = bit_interleave_zero(val0);
    val1 = bit_interleave_zero(val1);
    return val0 | (val1 << 1);
}

/// Bayer index for ordered dithering.
///
/// Let \f$ \mathbf{M} \in \mathbb{R}^{n,n} \f$ denote the Bayer index
/// matrix for dimension \f$ n = 2^B \f$, where
/// \f[
///      M_{i,j} =
///         \frac{1}{n^2}
///         \operatorname{reverse}_{2B}(
///         \operatorname{interleave}(\operatorname{xor}(j, i), i)).
/// \f]
///
/// This routine calculates \f$ n^2 M_{i,j} \f$, which is always an
/// integer.
///
template <std::size_t B, std::integral Int>
constexpr Int bayer_index(Int i, Int j) noexcept {
    static_assert(B > 0, "B must be greater than 0");
    static_assert(sizeof(Int) * 8 > B * 2, "Int does not have enough bits");
    i = i & ((1 << B) - 1);
    j = j & ((1 << B) - 1);
    return bit_reverse<Int>(bit_interleave<Int>(j ^ i, i)) >>
           (sizeof(Int) * 8 - B * 2);
}

/// Greatest common divisor.
template <std::integral Int>
constexpr Int gcd(Int a, Int b) noexcept {
    Int r = 0;
    while (a) {
        r = b % a;
        b = a;
        a = r;
    }
    return b;
}

/// Greatest common divisor with Bezout coefficients.
template <std::integral Int>
constexpr Int gcd_bezout(Int a, Int b, Int* x, Int* y) noexcept {
    Int s1 = 0, s0 = 1;
    Int t1 = 1, t0 = 0;
    Int r1 = b, r0 = a;
    while (r1 != 0) {
        Int q = r0 / r1;
        Int k;
        k = r1, r1 = r0 - q * r1, r0 = k;
        k = s1, s1 = s0 - q * s1, s0 = k;
        k = t1, t1 = t0 - q * t1, t0 = k;
    }
    if (x)
        *x = s0;
    if (y)
        *y = t0;
    return r0;
}

/// Least common multiple.
template <std::integral Int>
constexpr Int lcm(Int a, Int b) noexcept {
    if (!a || !b)
        return 0;
    else
        return a * b / pre::gcd(a, b);
}

/// Cantor pairing function.
///
/// \f[
///     \pi(x, y) =
///     \frac{1}{2}(x + y)(x + y + 1) + y
/// \f]
///
template <std::integral Int>
constexpr Int cantor(Int x, Int y) noexcept {
    if constexpr (!std::unsigned_integral<Int>) {
        // Natural numbers.
        std::make_unsigned_t<Int> xn = std::abs(x);
        std::make_unsigned_t<Int> yn = std::abs(y);
        xn = 2 * xn - (x > 0);
        yn = 2 * yn - (y > 0);
        return Int(cantor(xn, yn));
    }
    else {
        return ((x + y) * (x + y + 1)) / 2 + y;
    }
}

/// Cantor tuple pairing function.
///
/// \f[
///      \pi(x, y, z, \ldots) = \pi(\pi(x, y), z, \ldots)
/// \f]
///
template <std::integral Int, std::integral... Ints>
constexpr Int cantor(Int x, Int y, Int z, Ints... ws) noexcept {
    return cantor(cantor(x, y), z, ws...);
}

template <std::ranges::range Ints>
constexpr auto cantor(const Ints& ints) noexcept {
    using Int = std::ranges::range_value_t<Ints>;
    Int c = 0;
    for (Int x : ints)
        c = cantor(c, x);
    return c;
    static_assert(std::integral<Int>);
}

/** \} */

/// \name Floating point helpers
/** \{ */

/// Finite or alternative.
template <std::floating_point Float>
inline Float finite_or(Float x, Float x0) noexcept {
    return std::isfinite(x) ? x : x0;
}

/// Finite or zero.
template <std::floating_point Float>
inline Float finite_or_zero(Float x) noexcept {
    return finite_or(x, Float(0));
}

/// Wrap floating point number in range.
template <std::floating_point Float>
inline Float frepeat(
    Float x,
    std::type_identity_t<Float> a,
    std::type_identity_t<Float> b) noexcept {
    x -= a;
    b -= a;
    Float r = std::remainder(x, b);
    if (r < 0)
        r += b;
    return r + a;
}

/// Wrap floating point number in range and mirror with each repeat.
template <std::floating_point Float>
inline Float fmirror(
    Float x,
    std::type_identity_t<Float> a,
    std::type_identity_t<Float> b) noexcept {
    x -= a;
    b -= a;
    int q = 0;
    Float r = std::remquo(x, b, &q);
    if (r < 0) {
        r += b;
        q++;
    }
    if (q & 1)
        r = b - r;
    return r + a;
}

/// Increment float to next representable value.
template <std::floating_point Float>
inline Float float_incr(Float x) noexcept {
    if constexpr (
        std::same_as<Float, float> and //
        std::numeric_limits<float>::is_iec559) {
        std::uint32_t u = bit_cast<std::uint32_t>(x);
        if (u != (0x7f8UL << 20)) { // Not +Inf?
            if (u == (1UL << 31))   // Ignore -0.0.
                u = 0;
            if (u & (1UL << 31))
                u = u - 1;
            else
                u = u + 1;
        }
        return bit_cast<float>(u);
    }
    else if constexpr (
        std::same_as<Float, double> and
        std::numeric_limits<double>::is_iec559) {
        std::uint64_t u = bit_cast<std::uint64_t>(x);
        if (u != (0x7ffULL << 52)) { // Not +Inf?
            if (u == (1ULL << 63))   // Ignore -0.0.
                u = 0;
            if (u & (1ULL << 63))
                u = u - 1;
            else
                u = u + 1;
        }
        return bit_cast<double>(u);
    }
    else
        return std::nextafter(x, +Inf<Float>);
}

/// Decrement float to next representable value.
template <std::floating_point Float>
inline Float float_decr(Float x) noexcept {
    if constexpr (
        std::same_as<Float, float> and //
        std::numeric_limits<float>::is_iec559) {
        std::uint32_t u = bit_cast<std::uint32_t>(x);
        if (u != (0xff8UL << 20)) { // Not -Inf?
            if (u == 0)             // Ignore +0.0.
                u = (1UL << 31);
            if (u & (1UL << 31))
                u = u + 1;
            else
                u = u - 1;
        }
        return bit_cast<float>(u);
    }
    else if constexpr (
        std::same_as<Float, double> and
        std::numeric_limits<double>::is_iec559) {
        std::uint64_t u = bit_cast<std::uint64_t>(x);
        if (u != (0xfffULL << 52)) { // Not -Inf?
            if (u == 0)              // Ignore +0.0.
                u = (1ULL << 63);
            if (u & (1ULL << 63))
                u = u + 1;
            else
                u = u - 1;
        }
        return bit_cast<double>(u);
    }
    else
        return std::nextafter(x, -Inf<Float>);
}

/// Same signbit?
template <std::floating_point Float>
inline bool same_signbit(Float a, Float b) noexcept {
    return std::signbit(a) == std::signbit(b);
}

/// Is relatively tiny? (a much less than b)
template <std::floating_point Float>
inline bool is_tiny(Float a, Float b) noexcept {
    volatile Float x = a;
    volatile Float y = b;
    volatile Float z = x + y;
    return z == y;
}

/// Is relatively huge? (a much greater than b)
template <std::floating_point Float>
inline bool is_huge(Float a, Float b) noexcept {
    return is_tiny(b, a);
}

/// Quadratic roots.
template <std::floating_point Float>
inline bool quadratic(
    Float a, Float b, Float c, Float& x0, Float& x1) noexcept {
    if (is_tiny(a, b)) {
        x0 = x1 = -c / b;
        return std::isfinite(x0);
    }
    else {
        b /= a;
        c /= a;
        if (not std::isfinite(b) or //
            not std::isfinite(c))
            return false;
        Float d = b * b - 4 * c;
        if (not std::isfinite(d))
            d = b * (b - 4 * (c / b)); // Try again
        if (not std::isfinite(d) or d < 0)
            return false;
        x0 = -Float(0.5) * (b + std::copysign(std::sqrt(d), b));
        x1 = c / x0;
        if (x0 > x1)
            std::swap(x0, x1);
        return true;
    }
}

/// Polynomial value.
template <
    concepts::arithmetic_or_complex Coeff,
    concepts::arithmetic_or_complex Coord>
constexpr auto polyval(IteratorRange<const Coeff*> ps, Coord x) noexcept {
    std::decay_t<decltype(Coeff() * Coord())> y = {};
    for (Coeff p : ps)
        y = y * x + p;
    return y;
}

/// Rational polynomial value.
template <
    concepts::arithmetic_or_complex Coeff,
    concepts::arithmetic_or_complex Coord>
constexpr auto polyval(
    IteratorRange<const Coeff*> ps,
    IteratorRange<const Coeff*> qs,
    Coord x) noexcept {
    if constexpr (std::integral<Coord>)
        return polyval(ps, qs, double(x));
    else {
        ASSERT(ps.size() == qs.size());
        if (std::abs(x) <= 1)
            return polyval(ps, x) / polyval(qs, x);
        else {
            int n = ps.size();
            x = Coord(1) / x;
            std::decay_t<decltype(Coeff() * Coord())> yp = {};
            std::decay_t<decltype(Coeff() * Coord())> yq = {};
            for (int k = n - 1; k >= 0; k--) {
                yp = yp * x + ps[k];
                yq = yq * x + qs[k];
            }
            return yp / yq;
        }
    }
}

/// Fast inverse square root.
template <std::floating_point Float>
[[gnu::always_inline]] inline Float fast_inv_sqrt(Float x) noexcept {
    if constexpr (
        std::same_as<Float, float> and //
        pre::numeric_limits<float>::is_iec559) {
        float y = x;
        float h = y * 0.5f;
        std::uint32_t u = 0x5f375a86 - (bit_cast<std::uint32_t>(y) >> 1);
        y = bit_cast<float>(u);
        y = y * (1.5f - h * y * y);
        y = y * (1.5f - h * y * y);
        return y;
    }
    else {
        return 1 / pre::sqrt(x);
    }
}

/// Fast floor by int casting.
template <std::floating_point Float, std::integral Int = int>
constexpr Int fast_floor(Float x) noexcept {
    Int i(x);
    i = i - (Float(i) > x);
    return i;
}

/// Fast ceil by int casting.
template <std::floating_point Float, std::integral Int = int>
constexpr Int fast_ceil(Float x) noexcept {
    Int i(x);
    i = i + (Float(i) < x);
    return i;
}

/// Fast round by int casting.
template <std::floating_point Float, std::integral Int = int>
constexpr Int fast_round(Float x) noexcept {
    return fast_floor<Float, Int>(x + Float(0.5));
}

/// Fast trunc by int casting.
template <std::floating_point Float, std::integral Int = int>
constexpr Int fast_trunc(Float x) noexcept {
    return Int(x);
}

/// Packed cast.
///
/// This function converts between floating point values and (un)signed
/// normalized "packed" integer values.
///
template <concepts::arithmetic To, concepts::arithmetic From>
constexpr To packed_cast(From x) noexcept {
    if constexpr (std::floating_point<From>) {
        if constexpr (std::floating_point<To>)
            return x;
        else if constexpr (std::unsigned_integral<To>)
            return fast_round<From, To>(clamp(x, +0, +1) * Maximum<To>);
        else
            return fast_round<From, To>(clamp(x, -1, +1) * Maximum<To>);
    }
    else {
        if constexpr (std::same_as<To, From>)
            return x;
        else if constexpr (std::floating_point<To>)
            return To(x) / To(Maximum<From>);
        else
            return packed_cast<To>(packed_cast<double>(x));
    }
}

/// Fraction with respect to floor.
template <std::floating_point Float>
constexpr Float fract(Float x, int* i = nullptr) noexcept {
    int x0;
    if (i == nullptr)
        i = &x0;
    *i = fast_floor(x);
    return x - *i;
}

/// Calculate `sin` and `cos`, return as pair.
template <concepts::arithmetic Arith>
inline auto sincos(Arith x) noexcept {
    // On GCC, could manually call __builtin_sincos, but this should
    // be optimizable as is.
    return std::make_pair(std::sin(x), std::cos(x));
}

/// Multiply by \f$ \pi \f$, then take sine.
///
/// \note
/// This is equal to \f$ \sin(\pi x) \f$. However, this
/// exploits periodicity, using `remquo()` before multiplying by
/// \f$ \pi \f$ to be more accurate for large arguments.
///
template <std::floating_point Float>
inline Float sinpi(Float x) noexcept {
    int quo;
    Float rem = pre::remquo(x, Float(1), &quo);
    Float res = pre::sin(M_pi<Float> * rem);
    if (unsigned(quo) & 1)
        res = -res;
    return res;
}

/// Multiply by \f$ \pi \f$, then take cosine.
///
/// \note
/// This is equal to \f$ \cos(\pi x) \f$. However, this
/// exploits periodicity using `remquo()` before multiplying by
/// \f$ \pi \f$ to be more accurate for large arguments.
///
template <std::floating_point Float>
inline Float cospi(Float x) noexcept {
    int quo;
    Float rem = pre::remquo(x, Float(1), &quo);
    Float res = pre::cos(M_pi<Float> * rem);
    if (unsigned(quo) & 1)
        res = -res;
    return res;
}

/// Find minimum value of `sinpi()` on interval.
template <std::floating_point Float>
inline Float sinpi_min(Float x1, Float x2) noexcept {
    if (!(x1 < x2)) {
        std::swap(x1, x2);
    }
    int n1 = pre::fast_floor(2 * x1);
    int n2 = pre::fast_floor(2 * x2);
    if (n2 - n1 > 3) {
        return -1;
    }
    else {
        // 0 = min is -1
        // 1 = min is sinpi(x1)
        // 2 = min is sinpi(x2)
        // 3 = min is either sinpi(x1) or sinpi(x2)
        constexpr int lookup[4][4] = {
            {1, 3, 2, 0}, {0, 2, 2, 0}, {0, 0, 2, 0}, {1, 1, 3, 1}};
        n1 &= 3;
        n2 &= 3;
        switch (lookup[n1][n2]) {
        case 0: return -1;
        case 1: return pre::sinpi(x1);
        case 2: return pre::sinpi(x2);
        case 3: return pre::min(pre::sinpi(x1), pre::sinpi(x2));
        default: break;
        }
        // Unreachable.
        return 0;
    }
}

/// Find maximum value of `sinpi()` on interval.
template <std::floating_point Float>
inline Float sinpi_max(Float x1, Float x2) noexcept {
    if (!(x1 < x2)) {
        std::swap(x1, x2);
    }
    int n1 = pre::fast_floor(2 * x1);
    int n2 = pre::fast_floor(2 * x2);
    if (n2 - n1 > 3) {
        return +1;
    }
    else {
        // 0 = max is +1
        // 1 = max is sinpi(x1)
        // 2 = max is sinpi(x2)
        // 3 = max is either sinpi(x1) or sinpi(x2)
        constexpr int lookup[4][4] = {
            {2, 0, 0, 0}, {3, 1, 1, 1}, {2, 0, 1, 3}, {2, 0, 0, 2}};
        n1 &= 3;
        n2 &= 3;
        switch (lookup[n1][n2]) {
        case 0: return +1;
        case 1: return pre::sinpi(x1);
        case 2: return pre::sinpi(x2);
        case 3: return pre::max(pre::sinpi(x1), pre::sinpi(x2));
        default: break;
        }
        // Unreachable.
        return 0;
    }
}

/// Calculate `sinpi` and `cospi`, return as pair.
template <std::floating_point Float>
inline auto sincospi(Float x) noexcept {
    int quo;
    Float rem = pre::remquo(x, Float(1), &quo);
    auto [sinx, cosx] = pre::sincos(M_pi<Float> * rem);
    if (unsigned(quo) & 1) {
        sinx = -sinx;
        cosx = -cosx;
    }
    return std::make_pair(sinx, cosx);
}

/// Error function inverse.
template <std::floating_point Float>
inline Float erfinv(Float y) noexcept {
    Float w = -std::log((1 - y) * (1 + y));
    Float p;
    if (w < Float(5)) {
        w = w - Float(2.5);
        p = std::fma(w, Float(+2.81022636e-08), Float(+3.43273939e-7));
        p = std::fma(w, p, Float(-3.52338770e-6));
        p = std::fma(w, p, Float(-4.39150654e-6));
        p = std::fma(w, p, Float(+2.18580870e-4));
        p = std::fma(w, p, Float(-1.25372503e-3));
        p = std::fma(w, p, Float(-4.17768164e-3));
        p = std::fma(w, p, Float(+2.46640727e-1));
        p = std::fma(w, p, Float(+1.50140941));
    }
    else {
        w = std::sqrt(w) - 3;
        p = std::fma(w, Float(-2.00214257e-4), Float(+1.00950558e-4));
        p = std::fma(w, p, Float(+1.34934322e-3));
        p = std::fma(w, p, Float(-3.67342844e-3));
        p = std::fma(w, p, Float(+5.73950773e-3));
        p = std::fma(w, p, Float(-7.62246130e-3));
        p = std::fma(w, p, Float(+9.43887047e-3));
        p = std::fma(w, p, Float(+1.00167406));
        p = std::fma(w, p, Float(+2.83297682));
    }
    return p * y;
}

/** \} */

} // namespace pre
