/*-*- C++ -*-*/
#pragma once

namespace pre {

/// A static string.
///
/// \tparam Char
/// String character type.
///
/// \tparam BufSize
/// String character capacity, must be greater than 1 (because
/// 1 character is always used as null terminator).
///
template <typename Char, size_t BufSize>
struct StaticString : ArrayLike<StaticString<Char, BufSize>, Char> {
  public:
    // Sanity check.
    static_assert(std::integral<Char> && BufSize > 1);

    constexpr StaticString() noexcept = default;

    constexpr StaticString(const StaticString&) noexcept = default;

    constexpr StaticString(StaticString&&) noexcept = default;

    template <typename... Args>
    StaticString(const std::basic_string<Char, Args...>& str) {
        if (str.size() > max_size())
            throw std::invalid_argument(__func__);
        std::copy(str.begin(), str.end(), begin());
        len_ = str.size();
    }

    template <typename... Args>
    constexpr StaticString(const std::basic_string_view<Char, Args...>& str) {
        if (str.size() > max_size())
            throw std::invalid_argument(__func__);
        std::copy(str.begin(), str.end(), begin());
        len_ = str.size();
    }

    constexpr StaticString(const Char* str)
        : StaticString(std::basic_string_view<Char>(str)) {
    }

    template <size_t OtherBufSize>
    constexpr StaticString(const StaticString<Char, OtherBufSize>& str)
        : StaticString(str.view()) {
    }

    constexpr StaticString& operator=(const StaticString&) = default;

    constexpr StaticString& operator=(StaticString&&) = default;

  public:
    /// \name Container API
    /** \{ */

    constexpr size_t size() const noexcept {
        return len_;
    }

    constexpr size_t max_size() const noexcept {
        return BufSize - 1;
    }

    constexpr size_t capacity() const noexcept {
        return BufSize - 1;
    }

    constexpr Char* begin() noexcept {
        return &buf_[0];
    }

    constexpr const Char* begin() const noexcept {
        return &buf_[0];
    }

    constexpr Char* end() noexcept {
        return &buf_[0] + len_;
    }

    constexpr const Char* end() const noexcept {
        return &buf_[0] + len_;
    }

    constexpr void clear() noexcept {
        std::fill(begin(), end(), Char(0));
        len_ = 0;
    }

    constexpr void resize(size_t len) {
        if (len > max_size())
            throw std::invalid_argument(__func__);
        len_ = len;
        buf_[len] = Char(0);
    }

    /** \} */

    /// \name String
    /** \{ */

    constexpr size_t length() const noexcept {
        return len_;
    }

    constexpr const Char* c_str() const noexcept {
        return &buf_[0];
    }

    template <typename... Args>
    operator std::basic_string<Char, Args...>() const {
        return &buf_[0];
    }

    template <typename... Args>
    constexpr operator std::basic_string_view<Char, Args...>() const noexcept {
        return {&buf_[0], len_};
    }

    constexpr std::basic_string_view<Char> view() const noexcept {
        return *this;
    }

    template <size_t OtherBufSize>
    constexpr auto operator<=>(
            const StaticString<Char, OtherBufSize>& other) const noexcept {
        return view() <=> other.view();
    }

    template <size_t OtherBufSize>
    constexpr bool operator==(
            const StaticString<Char, OtherBufSize>& other) const noexcept {
        return view() == other.view();
    }

    /** \} */

  public:
    void serialize(auto& serializer) {
        serializer <=> buf_;
        serializer <=> len_;
    }

  private:
    Char buf_[BufSize] = {};

    size_t len_ = 0;
};

} // namespace pre
