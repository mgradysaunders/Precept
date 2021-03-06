/*-*- C++ -*-*/
#pragma once

namespace pre {

using ssize_t = std::make_signed_t<size_t>;

constexpr size_t operator""_KiB(unsigned long long x) noexcept {
    return x * 1024ULL;
}

constexpr size_t operator""_MiB(unsigned long long x) noexcept {
    return x * 1048576ULL;
}

constexpr size_t operator""_GiB(unsigned long long x) noexcept {
    return x * 1073741824ULL;
}

constexpr size_t operator""_KiB(long double x) noexcept {
    return x * 1024.0L;
}

constexpr size_t operator""_MiB(long double x) noexcept {
    return x * 1048576.0L;
}

constexpr size_t operator""_GiB(long double x) noexcept {
    return x * 1073741824.0L;
}

enum class Bytes : size_t {};

template <concepts::ostream Stream>
inline Stream& operator<<(Stream& stream, Bytes bytes) {
    auto prev_flags = stream.flags();
    auto prev_precision = stream.precision(2);
    stream.setf(std::ios_base::dec, std::ios_base::basefield);
    stream.setf(std::ios_base::fixed, std::ios_base::floatfield);
    stream.unsetf(std::ios_base::showpos);
    size_t sz = size_t(bytes);
    if (sz < 1_KiB)
        stream << sz << 'B';
    else if (sz < 1_MiB)
        stream << double(sz) / double(1_KiB) << "KiB";
    else if (sz < 1_GiB)
        stream << double(sz) / double(1_MiB) << "MiB";
    else
        stream << double(sz) / double(1_GiB) << "GiB";
    stream.flags(prev_flags);
    stream.precision(prev_precision);
    return stream;
}

template <size_t N>
struct sized_int {
    static constexpr auto select_type() noexcept {
        if constexpr (sizeof(std::int_least8_t) >= N)
            return std::int_least8_t();
        else if constexpr (sizeof(std::int_least16_t) >= N)
            return std::int_least16_t();
        else if constexpr (sizeof(std::int_least32_t) >= N)
            return std::int_least32_t();
        else if constexpr (sizeof(std::int_least64_t) >= N)
            return std::int_least64_t();
    }
    using type = std::decay_t<decltype(select_type())>;
};

template <size_t N>
struct sized_uint {
    static constexpr auto select_type() noexcept {
        if constexpr (sizeof(std::uint_least8_t) >= N)
            return std::uint_least8_t();
        else if constexpr (sizeof(std::uint_least16_t) >= N)
            return std::uint_least16_t();
        else if constexpr (sizeof(std::uint_least32_t) >= N)
            return std::uint_least32_t();
        else if constexpr (sizeof(std::uint_least64_t) >= N)
            return std::uint_least64_t();
    }
    using type = std::decay_t<decltype(select_type())>;
};

template <size_t N>
struct sized_float {
    static constexpr auto select_type() noexcept {
        if constexpr (sizeof(float) >= N)
            return 0.0f;
        else if constexpr (sizeof(double) >= N)
            return 0.0;
        else if constexpr (sizeof(long double) >= N)
            return 0.0L;
    }
    using type = std::decay_t<decltype(select_type())>;
};

template <size_t N>
using sized_int_t = typename sized_int<N>::type;

template <size_t N>
using sized_uint_t = typename sized_uint<N>::type;

template <size_t N>
using sized_float_t = typename sized_float<N>::type;

template <typename To, typename From>
inline To bit_cast(const From& from) noexcept
        requires(sizeof(To) == sizeof(From)) {
    To to;
    std::memcpy(&to, &from, sizeof(from));
    return to;
}

} // namespace pre
