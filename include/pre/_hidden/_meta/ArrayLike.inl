/*-*- C++ -*-*/
#pragma once

namespace pre {

/// Array-like CRTP.
template <typename Subclass, typename Value, bool NegativeFromBack = true>
struct ArrayLike {

    typedef size_t size_type;

    typedef std::ptrdiff_t difference_type;

    typedef Value value_type;

    typedef value_type& reference;

    typedef value_type* pointer;

    typedef const value_type& const_reference;

    typedef const value_type* const_pointer;

    typedef pointer iterator;

    typedef const_pointer const_iterator;

    typedef std::reverse_iterator<iterator> reverse_iterator;

    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    constexpr bool empty() const noexcept {
        return static_cast<const Subclass&>(*this).size() == 0;
    }

    constexpr const_iterator cbegin() const noexcept {
        return static_cast<const Subclass&>(*this).begin();
    }

    constexpr const_iterator cend() const noexcept {
        return static_cast<const Subclass&>(*this).end();
    }

    constexpr reverse_iterator rbegin() noexcept {
        return reverse_iterator(static_cast<Subclass&>(*this).end());
    }

    constexpr const_reverse_iterator rbegin() const noexcept {
        return const_reverse_iterator(
                static_cast<const Subclass&>(*this).end());
    }

    constexpr const_reverse_iterator crbegin() const noexcept {
        return rbegin(); // Force const
    }

    constexpr reverse_iterator rend() noexcept {
        return reverse_iterator(static_cast<Subclass&>(*this).begin());
    }

    constexpr const_reverse_iterator rend() const noexcept {
        return const_reverse_iterator(
                static_cast<const Subclass&>(*this).begin());
    }

    constexpr const_reverse_iterator crend() const noexcept {
        return rend(); // Force const
    }

    constexpr reference front() noexcept {
        return *static_cast<Subclass&>(*this).begin();
    }

    constexpr const_reference front() const noexcept {
        return *static_cast<const Subclass&>(*this).begin();
    }

    constexpr reference back() noexcept {
        return *rbegin();
    }

    constexpr const_reference back() const noexcept {
        return *rbegin();
    }

    constexpr pointer data() noexcept {
        return static_cast<Subclass&>(*this).begin();
    }

    constexpr const_pointer data() const noexcept {
        return static_cast<const Subclass&>(*this).begin();
    }

    template <std::integral Int>
    constexpr reference operator[](Int pos) noexcept {
        if constexpr (std::signed_integral<Int> && NegativeFromBack)
            if (pos < 0)
                pos += static_cast<Subclass&>(*this).size();
        return *(static_cast<Subclass&>(*this).begin() + pos);
    }

    template <std::integral Int>
    constexpr const_reference operator[](Int pos) const noexcept {
        if constexpr (std::signed_integral<Int> && NegativeFromBack)
            if (pos < 0)
                pos += static_cast<const Subclass&>(*this).size();
        return *(static_cast<const Subclass&>(*this).begin() + pos);
    }

    template <std::integral Int>
    constexpr reference at(Int pos) {
        Int sz = static_cast<Subclass&>(*this).size();
        if constexpr (std::signed_integral<Int> && NegativeFromBack)
            if (pos < 0)
                pos += sz;
        if (pos < 0 || pos >= sz)
            throw std::out_of_range(__func__);
        return *(static_cast<Subclass&>(*this).begin() + pos);
    }

    template <std::integral Int>
    constexpr const_reference at(Int pos) const {
        Int sz = static_cast<const Subclass&>(*this).size();
        if constexpr (std::signed_integral<Int> && NegativeFromBack)
            if (pos < 0)
                pos += sz;
        if (pos < 0 || pos >= sz)
            throw std::out_of_range(__func__);
        return *(static_cast<const Subclass&>(*this).begin() + pos);
    }
};

} // namespace pre
