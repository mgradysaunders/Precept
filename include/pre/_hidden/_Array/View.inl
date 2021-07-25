/*-*- C++ -*-*/
#pragma once

namespace pre {

/// Slice for array views.
struct Slice {
  public:
    constexpr Slice() noexcept = default;

    constexpr Slice(ssize_t a, ssize_t b) noexcept : from(a), to(b) {
    }

    constexpr Slice(ssize_t k) noexcept : from(k), to(k + 1) {
    }

    /// Index to slice from.
    ssize_t from = 0;

    /// Index to slice to (non-inclusive).
    ssize_t to = std::numeric_limits<ssize_t>::max();

  private:
    constexpr void do_slice(
            auto*& first, ssize_t& size, ssize_t& skip) noexcept {
        if (from < 0)
            from += size + 1;
        if (from > size)
            from = size;
        if (to < 0)
            to += size + 1;
        if (to > size)
            to = size;
        ssize_t diff = to - from;
        if (diff < 0) {
            --from;
            --to;
        }
        first += skip * from;
        size = diff < 0 ? -diff : diff;
        skip = diff < 0 ? -skip : skip;
    }

    template <typename, size_t>
    friend struct ArrayView;
};

namespace concepts {

template <typename T>
concept integral_or_slice =
        std::integral<std::decay_t<T>> || std::same_as<std::decay_t<T>, Slice>;

} // namespace concepts

enum SliceInt : ssize_t {};

constexpr Slice operator|(SliceInt a, SliceInt b) noexcept {
    return {ssize_t(a), ssize_t(b)};
}

template <std::integral Int>
constexpr Slice operator|(SliceInt a, Int b) noexcept {
    return {ssize_t(a), ssize_t(b)};
}

template <std::integral Int>
constexpr Slice operator|(Int a, SliceInt b) noexcept {
    return {ssize_t(a), ssize_t(b)};
}

template <typename Value, size_t Rank>
struct ArrayView_initializer_list {
    using type = std::initializer_list<
            typename ArrayView_initializer_list<Value, Rank - 1>::type>;
};

template <typename Value>
struct ArrayView_initializer_list<Value, 0> {
    using type = Value;
};

/// An array view.
template <typename Value, size_t Rank>
struct ArrayView {
  public:
    static_assert(Rank > 0);

    struct Iterator {
      public:
        typedef std::ptrdiff_t difference_type;

        typedef std::
                conditional_t<Rank == 1, Value, ArrayView<Value, Rank - 1>>
                        value_type;

        typedef std::
                conditional_t<Rank == 1, Value&, ArrayView<Value, Rank - 1>>
                        reference;

        typedef std::conditional_t<Rank == 1, Value*, void> pointer;

        typedef std::random_access_iterator_tag iterator_category;

      private:
        constexpr Iterator() noexcept = default;

        [[gnu::nonnull]] constexpr Iterator(
                Value* first,
                const ssize_t* sizes,
                const ssize_t* skips) noexcept
            : first_(first), sizes_(sizes), skips_(skips) {
        }

      public:
        constexpr Iterator& operator++() noexcept {
            first_ += *skips_;
            return *this;
        }

        constexpr Iterator& operator--() noexcept {
            first_ -= *skips_;
            return *this;
        }

        constexpr Iterator operator++(int) noexcept {
            const Iterator copy = *this;
            operator++();
            return copy;
        }

        constexpr Iterator operator--(int) noexcept {
            const Iterator copy = *this;
            operator--();
            return copy;
        }

        constexpr Iterator& operator+=(difference_type count) noexcept {
            first_ += *skips_ * count;
            return *this;
        }

        constexpr Iterator operator+(difference_type count) const noexcept {
            return Iterator(*this) += count;
        }

        friend constexpr Iterator operator+(
                difference_type count, Iterator itr) noexcept {
            return itr + count;
        }

        constexpr Iterator& operator-=(difference_type count) noexcept {
            first_ -= *skips_ * count;
            return *this;
        }

        constexpr Iterator operator-(difference_type count) const noexcept {
            return Iterator(*this) -= count;
        }

        constexpr difference_type operator-(Iterator other) const noexcept {
            return (first_ - other.first_) / *other.skips_;
        }

        constexpr auto operator<=>(const Iterator& other) const noexcept {
            if (*skips_ > 0) {
                if (first_ < other.first_)
                    return std::strong_ordering::less;
                if (first_ > other.first_)
                    return std::strong_ordering::greater;
            }
            else {
                if (first_ > other.first_)
                    return std::strong_ordering::less;
                if (first_ < other.first_)
                    return std::strong_ordering::greater;
            }
            return std::strong_ordering::equal;
        }

        constexpr bool operator==(const Iterator&) const noexcept = default;

        constexpr reference operator*() noexcept {
            if constexpr (Rank == 1)
                return *first_;
            else
                return {first_, sizes_ + 1, skips_ + 1};
        }

      private:
        Value* first_ = nullptr;
        const ssize_t* sizes_ = nullptr;
        const ssize_t* skips_ = nullptr;

        friend struct ArrayView;
    };

    using iterator = Iterator;

    using initializer_list =
            typename ArrayView_initializer_list<Value, Rank>::type;

  public:
    constexpr ArrayView() noexcept = default;

    template <concepts::random_access_range_of<Value> Range>
    constexpr ArrayView(Range& range) requires(Rank == 1) {
        first = &*std::ranges::begin(range);
        sizes[0] = std::ranges::size(range);
        skips[0] = 1;
    }

    constexpr ArrayView(Value* ptr, ssize_t size0, ssize_t skip0 = 1) noexcept
            requires(Rank == 1) {
        first = ptr;
        sizes[0] = size0;
        skips[0] = skip0;
    }

    [[gnu::nonnull]] constexpr ArrayView(
            Value* ptr, const ssize_t* psizes, const ssize_t* pskips) noexcept
        : first(ptr) {
        std::copy(psizes, psizes + Rank, &sizes[0]);
        std::copy(pskips, pskips + Rank, &skips[0]);
    }

    constexpr ArrayView(Value* ptr, const ArrayIndex<Rank>& sz) noexcept
        : first(ptr), sizes(sz) {
        // Default implementation.
        skips[Rank - 1] = 1;
        for (size_t dim = Rank - 1; dim > 0; --dim)
            skips[dim - 1] = skips[dim] * sizes[dim];
    }

    constexpr ArrayView(
            Value* ptr,
            const ArrayIndex<Rank>& sz,
            const ArrayIndex<Rank>& sk) noexcept
        : first(ptr), sizes(sz), skips(sk) {
    }

    constexpr ArrayView(const ArrayView&) noexcept = default;

    constexpr ArrayView(ArrayView&&) noexcept = default;

    constexpr ArrayView& operator=(const ArrayView&) noexcept = default;

    constexpr ArrayView& operator=(ArrayView&&) noexcept = default;

    constexpr ArrayView& operator=(Value all) noexcept {
        for (auto itr = begin(); itr != end(); ++itr)
            *itr = all;
        return *this;
    }

    constexpr ArrayView& operator=(initializer_list list) noexcept {
        ASSERT(size() == ssize_t(list.size()));
        auto itr1 = list.begin();
        auto itr2 = begin();
        for (; itr2 != end(); ++itr1, ++itr2)
            *itr2 = *itr1;
        return *this;
    }

    template <typename Func>
    constexpr ArrayView& operator=(LazyArray<Func, Rank>&& lazy) noexcept {
        ASSERT(sizes == lazy.sizes);
        for (auto ind : ArrayRange(sizes))
            operator[](ind) = std::forward<LazyArray<Func, Rank>>(lazy)(ind);
        return *this;
    }

  public:
    constexpr auto rank() const noexcept {
        return Rank;
    }

    constexpr auto size() const noexcept {
        return sizes[0];
    }

    constexpr bool empty() const noexcept {
        return sizes.prod() == 0;
    }

    constexpr Iterator begin() noexcept {
        return {first, &sizes[0], &skips[0]};
    }

    constexpr Iterator end() noexcept {
        return begin() + sizes[0];
    }

  public:
    /// \name Indexing
    ///
    /// Views support indexing with the ordinary `operator[]` with an
    /// integer argument, such that the values can be accessed with the
    /// same syntax as built-in arrays, for example `v[0][1][2]`. It is
    /// also possible to use `operator()` instead like `v(0, 1, 2)`.
    ///
    /// _Slicing_ refers to the process of forming a range or subset
    /// of a view in one or more arbitrary dimensions. Suppose `v`
    /// is a rank-1 view and we want to form a smaller rank-1 view that
    /// only contains values from index 2 (inclusive) up to index 5
    /// (non-inclusive), then what we are looking for is `v[Slice(2, 5)]`.
    /// If we want to slice column 1 of a rank-2 view `m` that is 4 rows
    /// by 3 columns, then we would use `operator()` with a slice in the
    /// first dimension `m(Slice(0, 4), 1)`. It is also possible to
    /// slice along more than one axis at once, to obtain a block in a
    /// matrix for example.
    ///
    /// \note
    /// It necessary to use `operator()` when slicing along an axis that
    /// is not the last axis. Notice that `m[Slice(0, 4)][1]` does _not_
    /// return column 1 because it is actually two applications of
    /// `operator[]`.
    ///
    /** \{ */

    template <size_t R>
    constexpr decltype(auto) operator[](ArrayIndex<R> k) noexcept {
        static_assert(R <= Rank);
        if constexpr (R == Rank)
            return (*(first + k.linearize(skips)));
        else {
            ArrayView<Value, Rank - R> result;
            result.first = first + k.linearize(skips);
            std::copy(sizes.begin() + R, sizes.end(), result.sizes.begin());
            std::copy(skips.begin() + R, skips.end(), result.skips.begin());
            return result;
        }
    }

    template <size_t R>
    constexpr decltype(auto) operator[](ArrayIndex<R> k) const noexcept {
        static_assert(R <= Rank);
        if constexpr (R == Rank)
            return (*(first + k.linearize(skips)));
        else {
            ArrayView<Value, Rank - R> result;
            result.first = first + k.linearize(skips);
            std::copy(sizes.begin() + R, sizes.end(), result.sizes.begin());
            std::copy(skips.begin() + R, skips.end(), result.skips.begin());
            return result;
        }
    }

    template <concepts::integral_or_slice P>
    constexpr decltype(auto) operator[](P p) noexcept {
        if constexpr (std::integral<std::decay_t<P>>) {
            if constexpr (std::signed_integral<std::decay_t<P>>)
                if (p < 0)
                    p += sizes[0];
            return *(begin() + p);
        }
        else {
            ArrayView<Value, Rank> result = *this;
            Slice(p).do_slice(
                    result.first,    //
                    result.sizes[0], //
                    result.skips[0]);
            return result;
        }
    }

    template <concepts::integral_or_slice P, concepts::integral_or_slice... Q>
    constexpr decltype(auto) operator()(P p, Q&&... q) noexcept {
        constexpr size_t IndexCount = 1 + sizeof...(Q);
        constexpr size_t SliceCount =
                (size_t(std::same_as<std::decay_t<P>, Slice>) + ... +
                 size_t(std::same_as<std::decay_t<Q>, Slice>));
        if constexpr (SliceCount == 0) {
            if constexpr (sizeof...(Q) == 0)
                return operator[](p);
            else
                return operator[](p).operator()(std::forward<Q>(q)...);
        }
        else {
            Slice slices[] = {p, std::forward<Q>(q)...};
            constexpr bool is_slice[] = {
                    not std::integral<std::decay_t<P>>,
                    not std::integral<std::decay_t<Q>>...};
            ArrayView<Value, Rank - IndexCount + SliceCount> result;
            result.first = first;
            for (size_t dim = 0, off = 0; dim < IndexCount; dim++) {
                if (is_slice[dim]) {
                    result.sizes[off] = sizes[dim];
                    result.skips[off] = skips[dim];
                    slices[dim].do_slice(
                            result.first,      //
                            result.sizes[off], //
                            result.skips[off]);
                    off++;
                }
                else {
                    if (slices[dim].from < 0)
                        slices[dim].from += sizes[dim];
                    result.first += skips[dim] * slices[dim].from;
                }
            }
            std::copy(
                    sizes.begin() + IndexCount, //
                    sizes.end(), result.sizes.begin() + SliceCount);
            std::copy(
                    skips.begin() + IndexCount, //
                    skips.end(), result.skips.begin() + SliceCount);
            return result;
        }
    }

    /** \} */

  public:
    /// \name Lazyification
    /** \{ */
    // clang-format off

    [[gnu::always_inline]] 
    constexpr auto lazy() noexcept {
        return LazyArray(
                [=, *this](auto ind) constexpr noexcept {
                    return operator[](ind);
                },
                sizes);
    }

    [[gnu::always_inline]]
    constexpr auto operator*() noexcept {
        return lazy();
    }

    // clang-format on
    /** \} */

  public:
    /// \name Utilities
    /** \{ */

    template <typename Func>
    constexpr void for_each(Func&& func) noexcept {
        for (auto ind : ArrayRange(sizes))
            std::invoke(std::forward<Func>(func), operator[](ind));
    }

    constexpr void get_each(std::output_iterator<Value> auto itr) noexcept {
        for (auto ind : ArrayRange(sizes))
            *itr++ = operator[](ind);
    }

    constexpr void set_each(std::input_iterator auto itr) noexcept {
        for (auto ind : ArrayRange(sizes))
            operator[](ind) = *itr++;
    }

    template <typename Other>
    constexpr void swap_each(ArrayView<Other, Rank> other) noexcept {
        auto& lhs = *this;
        auto& rhs = other;
        ASSERT(lhs.sizes == rhs.sizes);
        for (auto ind : ArrayRange(sizes))
            std::swap(lhs[ind], rhs[ind]);
    }

    /// Any bottom-level values true?
    constexpr bool any() noexcept {
        for (auto each : *this)
            if constexpr (Rank == 1) {
                if (each)
                    return true;
            }
            else {
                if (each.any())
                    return true;
            }
        return false;
    }

    /// All bottom-level values true?
    constexpr bool all() noexcept {
        for (auto each : *this)
            if constexpr (Rank == 1) {
                if (not each)
                    return false;
            }
            else {
                if (not each.all())
                    return false;
            }
        return true;
    }

    /// Reorder.
    template <std::integral... Ints>
    constexpr ArrayView reorder(Ints... ints) const noexcept {
        static_assert(sizeof...(Ints) == Rank);
        ArrayView result{
                first,                  //
                sizes.swizzle(ints...), //
                skips.swizzle(ints...)};
        ASSERT(std::is_permutation(
                result.skips.begin(), //
                result.skips.end(), skips.begin()));
        return result;
    }

    /** \} */

  public:
    /// \name Utilities (Rank-1)
    /** \{ */

    constexpr ssize_t argmin() noexcept requires(Rank == 1) {
        return std::min_element(begin(), end()) - begin();
    }

    constexpr ssize_t argmax() noexcept requires(Rank == 1) {
        return std::max_element(begin(), end()) - begin();
    }

    constexpr auto sum() noexcept requires(Rank == 1) {
        std::decay_t<Value> res = 0;
        for (Value each : *this)
            res += each;
        return res;
    }

    constexpr auto prod() noexcept requires(Rank == 1) {
        std::decay_t<Value> res = 1;
        for (Value each : *this)
            res *= each;
        return res;
    }

    /** \} */

  public:
    /// \name Utilities (Rank-2)
    /** \{ */

    constexpr bool is_square() const noexcept requires(Rank == 2) {
        return sizes[0] == sizes[1];
    }

    constexpr ssize_t rows() noexcept requires(Rank == 2) {
        return sizes[0];
    }

    constexpr ssize_t cols() noexcept requires(Rank == 2) {
        return sizes[1];
    }

    constexpr ArrayView<Value, 1> row(ssize_t p) noexcept requires(Rank == 2) {
        return (*this)[p];
    }

    constexpr ArrayView<Value, 1> col(ssize_t p) noexcept requires(Rank == 2) {
        return (*this)(Slice(), p);
    }

    constexpr void swap_rows(ssize_t k0, ssize_t k1) noexcept
            requires(Rank == 2) {
        if (k0 != k1)
            row(k0).swap_each(row(k1));
    }

    constexpr void swap_cols(ssize_t k0, ssize_t k1) noexcept
            requires(Rank == 2) {
        if (k0 != k1)
            col(k0).swap_each(col(k1));
    }

    constexpr ArrayView<Value, 1> diag(ssize_t p = 0) noexcept
            requires(Rank == 2) {
        ssize_t size0 = sizes[0], size1 = sizes[1];
        ssize_t skip0 = skips[0], skip1 = skips[1];
        if (!(first && p > -size0 && p < +size1))
            return {};
        if (p < 0) {
            p = -p;
            std::swap(size0, size1);
            std::swap(skip0, skip1);
        }
        return {first + skip1 * p, std::min(size0, size1 - p), skip0 + skip1};
    }

    constexpr auto trace() noexcept requires(Rank == 2) {
        return diag().sum();
    }

    constexpr ArrayView transpose() noexcept requires(Rank == 2) {
        ArrayView result = *this;
        std::swap(result.sizes[0], result.sizes[1]);
        std::swap(result.skips[0], result.skips[1]);
        return result;
    }

    constexpr ArrayView transpose(ssize_t k0, ssize_t k1) noexcept
            requires(Rank > 2) {
        ArrayView result = *this;
        std::swap(result.sizes[k0], result.sizes[k1]);
        std::swap(result.skips[k0], result.skips[k1]);
        return result;
    }

    /** \} */

  public:
    /// Implicit cast as bool. (Not empty?)
    constexpr operator bool() const noexcept {
        return !empty();
    }

    /// Implicit cast as const.
    constexpr operator ArrayView<const Value, Rank>() const noexcept {
        return {first, &sizes[0], &skips[0]};
    }

  public:
    Value* first = nullptr;
    ArrayIndex<Rank> sizes = {};
    ArrayIndex<Rank> skips = {};

  public:
    /// Read from `std::basic_istream`.
    ///
    /// Format is `[v0,v1,...]`.
    /// Sets `std::ios_base::failbit` on error.
    ///
    template <concepts::istream Stream>
    friend Stream& operator>>(Stream& stream, ArrayView arr) {
        using Char = typename Stream::char_type;
        using CharTraits = typename Stream::traits_type;
        auto consume = [&](char what) {
            Char ch;
            if (!(stream >> ch) ||
                !CharTraits::eq(ch, CharTraits::to_char_type(what))) {
                stream.setstate(std::ios_base::failbit);
                return false;
            }
            return true;
        };
        if (consume('[')) {
            if (!arr.empty()) {
                auto itr = arr.begin();
                stream >> *itr++;
                while (itr != arr.end()) {
                    if (!consume(','))
                        return stream;
                    stream >> *itr++;
                }
            }
            consume(']');
        }
        return stream;
    }

    /// Write into `std::basic_ostream`.
    ///
    /// Format is `[v0,v1,...]`.
    ///
    template <concepts::ostream Stream>
    friend Stream& operator<<(Stream& stream, ArrayView arr) {
        stream << '[';
        if (!arr.empty()) {
            auto itr = arr.begin();
            stream << *itr++;
            while (itr != arr.end()) {
                stream << ',';
                stream << *itr++;
            }
        }
        stream << ']';
        return stream;
    }
};

template <typename Value, size_t Rank>
ArrayView(Value*, ArrayIndex<Rank>) -> ArrayView<Value, Rank>;

template <typename Value, std::integral Int>
ArrayView(Value*, Int) -> ArrayView<Value, 1>;

template <std::ranges::random_access_range Range>
ArrayView(Range&) -> ArrayView<std::ranges::range_value_t<Range>, 1>;

template <std::ranges::random_access_range Range>
ArrayView(const Range&)
        -> ArrayView<const std::ranges::range_value_t<Range>, 1>;

} // namespace pre

#include "View_operators.inl"
