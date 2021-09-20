/*-*- C++ -*-*/
#pragma once

namespace pre {

/// A multi-dimensional array index helper.
template <size_t Rank>
struct ArrayIndex : ArrayLike<ArrayIndex<Rank>, ssize_t> {
  public:
    constexpr ArrayIndex() noexcept = default;

    template <std::integral... Ints>
    constexpr ArrayIndex(Ints... ints) noexcept : values{ssize_t(ints)...} {
        static_assert(sizeof...(Ints) == Rank);
    }

    template <std::ranges::input_range Ints>
    constexpr ArrayIndex(const Ints& arr) noexcept {
        std::copy_n(
            std::ranges::begin(arr),
            std::min(Rank, size_t(std::ranges::size(arr))), &values[0]);
    }

    template <std::input_iterator Iterator, size_t IgnoreRank>
    constexpr ArrayIndex(
        Iterator& itr, const ArrayIndex<IgnoreRank>& ignore) noexcept {
        for (size_t dim = 0; dim < Rank; dim++)
            if (not ignore.contains(dim))
                values[dim] = *itr++;
    }

  public:
    /// \name Container API
    /** \{ */

    constexpr size_t size() const noexcept {
        return Rank;
    }

    constexpr ssize_t* begin() noexcept {
        return &values[0];
    }

    constexpr const ssize_t* begin() const noexcept {
        return &values[0];
    }

    constexpr ssize_t* end() noexcept {
        return &values[0] + Rank;
    }

    constexpr const ssize_t* end() const noexcept {
        return &values[0] + Rank;
    }

    constexpr bool contains(ssize_t value) const noexcept {
        return std::find(begin(), end(), value) != end();
    }

    /** \} */

  public:
    /// Swizzle/permute.
    template <std::integral... Ints>
    constexpr auto swizzle(Ints... ints) const noexcept {
        return ArrayIndex<sizeof...(Ints)>{this->operator[](ints)...};
    }

    /// Increment, wrapping according to sizes.
    constexpr void increment(const ArrayIndex& sizes) noexcept {
        auto itr1 = sizes.rbegin();
        auto itr2 = this->rbegin();
        for (; itr2 < this->rend(); ++itr1, ++itr2) {
            ++*itr2;
            if (*itr2 >= *itr1)
                *itr2 = 0;
            else
                return;
        }
    }

    /// Decrement, wrapping according to sizes.
    constexpr void decrement(const ArrayIndex& sizes) noexcept {
        auto itr1 = sizes.rbegin();
        auto itr2 = this->rbegin();
        for (; itr2 < this->rend(); ++itr1, ++itr2) {
            --*itr2;
            if (*itr2 < 0)
                *itr2 = *itr1 - 1;
            else
                return;
        }
    }

    /// Linearize.
    constexpr ssize_t linearize(const ArrayIndex& skips) const noexcept {
        ssize_t res = 0;
        auto itr1 = skips.begin();
        auto itr2 = this->begin();
        for (; itr2 < this->end(); ++itr1, ++itr2)
            res += (*itr1) * (*itr2);
        return res;
    }

    /// Compute product.
    constexpr ssize_t prod() const noexcept {
        ssize_t res = values[0];
        for (size_t dim = 1; dim < Rank; dim++)
            res *= values[dim];
        return res;
    }

    template <size_t HeadRank>
    constexpr ArrayIndex<HeadRank> head() const noexcept {
        return IteratorRange(begin(), begin() + HeadRank);
        static_assert(HeadRank <= Rank);
    }

    template <size_t TailRank>
    constexpr ArrayIndex<TailRank> tail() const noexcept {
        return IteratorRange(begin() + Rank - TailRank, end());
        static_assert(TailRank <= Rank);
    }

    /// Join/concatenate with other ints.
    template <std::integral... Ints>
    constexpr auto join(Ints... other) const noexcept {
        return join(ArrayIndex<sizeof...(Ints)>(other...));
    }

    /// Join/concatenate with other.
    template <size_t JoinRank>
    constexpr auto join(const ArrayIndex<JoinRank>& other) const noexcept {
        ArrayIndex<Rank + JoinRank> res;
        std::copy(this->begin(), this->end(), res.begin());
        std::copy(other.begin(), other.end(), res.begin() + Rank);
        return res;
    }

    /// Omit dimensions.
    template <size_t OmitRank>
    constexpr auto omit(const ArrayIndex<OmitRank>& which) const noexcept {
        ArrayIndex<Rank - OmitRank> kept;
        ArrayIndex<Rank> omitted;
        auto itr1 = kept.begin();
        auto itr2 = omitted.begin();
        for (size_t dim = 0; dim < Rank; dim++) {
            if (not which.contains(dim))
                *itr1++ = values[dim];
            else
                *itr2++ = values[dim];
        }
        return std::make_pair(kept, omitted);
        static_assert(Rank >= OmitRank);
    }

    /// Levi-Civita permutation sign.
    constexpr int levi_civita() const noexcept {
        if constexpr (Rank == 0)
            return 0;
        else if constexpr (Rank == 1)
            return values[0] == 0;
        else if constexpr (Rank == 2) {
            if (*this == ArrayIndex{0, 1})
                return +1;
            if (*this == ArrayIndex{1, 0})
                return -1;
            return 0;
        }
        else {
            int visit[Rank] = {};
            for (ssize_t value : values)
                if (0 <= value and value < ssize_t(Rank))
                    visit[value]++;
            for (ssize_t value : values)
                if (visit[value] != 1)
                    return 0;
            std::fill(visit, visit + Rank, 0);
            int count = 0; // Count even cycles
            for (size_t dim = 0; dim < Rank; dim++)
                if (visit[dim] == 0) {
                    ssize_t len = 0;
                    ssize_t pos = dim;
                    for (; visit[pos] == 0; pos = values[pos], len++)
                        visit[pos] = 1;
                    count += ((len & 1) == 0);
                }
            return (count & 1) ? -1 : +1;
        }
    }

    constexpr auto operator<=>(const ArrayIndex& other) const noexcept {
        auto itr1 = this->begin();
        auto itr2 = other.begin();
        for (; itr2 < other.end(); ++itr1, ++itr2) {
            if (*itr1 < *itr2)
                return std::strong_ordering::less;
            if (*itr1 > *itr2)
                return std::strong_ordering::greater;
        }
        return std::strong_ordering::equal;
    }

    constexpr bool operator==(const ArrayIndex& other) const noexcept {
        auto itr1 = this->begin();
        auto itr2 = other.begin();
        for (; itr2 < other.end(); ++itr1, ++itr2)
            if (*itr1 != *itr2)
                return false;
        return true;
    }

  public:
    ssize_t values[Rank > 1 ? Rank : 1] = {};
};

// template <>
// struct ArrayIndex<0> {};

template <std::integral... Ints>
ArrayIndex(Ints... ints) -> ArrayIndex<sizeof...(Ints)>;

template <std::integral Int, size_t Rank>
ArrayIndex(const std::array<Int, Rank>&) -> ArrayIndex<Rank>;

template <std::integral Int, size_t Rank>
ArrayIndex(const Array<Int, Rank>&) -> ArrayIndex<Rank>;

template <size_t Rank>
struct ArrayRange {
    struct IteratorSentinel {};

    struct Iterator {
        typedef ArrayIndex<Rank> value_type;

        typedef ArrayIndex<Rank> reference;

        typedef std::input_iterator_tag iterator_category;

        constexpr Iterator() noexcept = default;

        constexpr Iterator(ArrayIndex<Rank> lim) noexcept
            : count(lim.prod()), limit(lim) {
        }

        constexpr Iterator& operator++() noexcept {
            index.increment(limit);
            count--;
            return *this;
        }

        constexpr ArrayIndex<Rank> operator*() const noexcept {
            return index;
        }

        constexpr bool operator==(IteratorSentinel) const noexcept {
            return count <= 0;
        }

        constexpr bool operator!=(IteratorSentinel) const noexcept {
            return count > 0;
        }

        ssize_t count = 0;

        ArrayIndex<Rank> limit = {};

        ArrayIndex<Rank> index = {};
    };

    constexpr ArrayRange() noexcept = default;

    template <typename... Args>
    constexpr ArrayRange(Args&&... args) noexcept
        : limit(std::forward<Args>(args)...) {
    }

    constexpr Iterator begin() const noexcept {
        return {limit};
    }

    constexpr IteratorSentinel end() const noexcept {
        return {};
    }

    ArrayIndex<Rank> limit;
};

template <size_t Rank>
ArrayRange(ArrayIndex<Rank>) -> ArrayRange<Rank>;

template <std::integral... Ints>
ArrayRange(Ints... ints) -> ArrayRange<sizeof...(Ints)>;

template <std::integral Int, size_t Rank>
ArrayRange(const std::array<Int, Rank>&) -> ArrayRange<Rank>;

template <std::integral Int, size_t Rank>
ArrayRange(const Array<Int, Rank>&) -> ArrayRange<Rank>;

} // namespace pre
