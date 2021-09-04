/*-*- C++ -*-*/
#pragma once

namespace pre {

/// An iterator range.
template <std::input_or_output_iterator Iterator>
class IteratorRange {
  public:
    typedef typename std::iterator_traits<Iterator>::value_type value_type;

    typedef typename std::iterator_traits<Iterator>::pointer pointer;

    typedef typename std::iterator_traits<Iterator>::reference reference;

    typedef typename std::iterator_traits<Iterator>::difference_type
        difference_type;

  public:
    constexpr IteratorRange() = default;

    constexpr IteratorRange(Iterator from, Iterator to)
        : begin_(from), end_(to) {
    }

    constexpr IteratorRange(const std::pair<Iterator, Iterator>& from_to)
        : begin_(from_to.first), end_(from_to.second) {
    }

    template <std::ranges::range Range>
    constexpr IteratorRange(Range&& range)
        : begin_(std::ranges::begin(range)), end_(std::ranges::end(range)) {
    }

    constexpr IteratorRange(const IteratorRange&) = default;

    constexpr IteratorRange(IteratorRange&&) = default;

    constexpr IteratorRange& operator=(const IteratorRange&) = default;

    constexpr IteratorRange& operator=(IteratorRange&&) = default;

    constexpr Iterator begin() const {
        return begin_;
    }

    constexpr Iterator end() const {
        return end_;
    }

    constexpr difference_type size() const {
        return empty() ? 0 : std::distance(begin_, end_);
    }

    constexpr bool empty() const {
        return begin_ == end_;
    }

    constexpr reference front() {
        return operator[](0);
    }

    constexpr reference back() {
        return operator[](size() - 1);
    }

    constexpr reference operator[](difference_type pos) {
        return *std::next(begin_, pos);
    }

    template <typename Value>
    constexpr void fill(const Value& value) {
        std::fill(begin_, end_, value);
    }

    template <typename Value>
    constexpr bool contains(const Value& value) {
        return std::find(begin_, end_, value) != end_;
    }

    template <typename Other>
    constexpr operator IteratorRange<Other>() {
        return {Other(begin_), Other(end_)};
    }

    template <typename Other>
    constexpr auto operator|(Other&& other) const noexcept {
        return std::ranges::subrange(begin_, end_) |
               std::forward<Other>(other);
    }

  private:
    Iterator begin_ = {};

    Iterator end_ = {};
};

template <std::input_or_output_iterator Iterator>
IteratorRange(const Iterator&, const Iterator&) -> IteratorRange<Iterator>;

template <std::input_or_output_iterator Iterator>
IteratorRange(const std::pair<Iterator, Iterator>&) -> IteratorRange<Iterator>;

template <std::ranges::range Range>
IteratorRange(Range&&) -> IteratorRange<std::ranges::iterator_t<Range>>;

} // namespace pre
