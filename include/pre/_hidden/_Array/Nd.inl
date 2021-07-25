/*-*- C++ -*-*/
#pragma once

namespace pre {

template <typename Value, typename Alloc = std::allocator<Value>>
class NdArray {
  public:
    /// \name Container typedefs
    /** \{ */

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

    /** \} */

  public:
    NdArray() = default;

    NdArray(const NdArray&) = default;

    NdArray(NdArray&&) = default;

    NdArray& operator=(const NdArray&) = default;

    NdArray& operator=(NdArray&&) = default;

    // TODO

  public:
    /// \name Container API
    /** \{ */

    void clear() {
        vals_.clear();
        dims_.clear();
    }

    template <std::ranges::input_range Dims>
    void resize(Dims&& dims) {
        dims_.resize(std::ranges::size(dims));
        std::copy(
                std::ranges::begin(dims),
                std::ranges::end(dims), //
                dims_.begin());
        ssize_t total_size = 1;
        for (ssize_t& dim : dims_) {
            if (dim < 0)
                dim = 0;
            total_size *= dim;
        }
        if (total_size == 0) {
            dims_.clear();
            vals_.clear();
        }
        else
            vals_.resize(total_size);
    }

    template <std::integral... Dims>
    void resize(Dims... dims) {
        this->resize(std::array{dims...});
    }

    template <typename... Args>
    void resize_like(const NdArray<Args...>& other) {
        this->resize(other.sizes());
    }

    size_t rank() const noexcept {
        return dims_.size();
    }

    size_t size() const noexcept {
        return vals_.size();
    }

    auto sizes() const noexcept {
        return IteratorRange(dims_.begin(), dims_.end());
    }

    template <size_t Rank>
    ArrayIndex<Rank> sizes() const {
        if (rank() > Rank)
            throw std::logic_error(__func__);
        ArrayIndex<Rank> dims;
        std::fill(dims.begin(), dims.end(), 1);
        std::copy(dims_.begin(), dims_.end(), dims.begin());
        return dims;
    }

    bool empty() const noexcept {
        return vals_.empty();
    }

    iterator begin() noexcept {
        return vals_.data();
    }

    const_iterator begin() const noexcept {
        return vals_.data();
    }

    const_iterator cbegin() const noexcept {
        return begin(); // Force const
    }

    iterator end() noexcept {
        return begin() + size();
    }

    const_iterator end() const noexcept {
        return begin() + size();
    }

    const_iterator cend() const noexcept {
        return end(); // Force const
    }

    reverse_iterator rbegin() noexcept {
        return reverse_iterator(end());
    }

    const_reverse_iterator rbegin() const noexcept {
        return const_reverse_iterator(end());
    }

    const_reverse_iterator crbegin() const noexcept {
        return rbegin(); // Force const
    }

    reverse_iterator rend() noexcept {
        return reverse_iterator(begin());
    }

    const_reverse_iterator rend() const noexcept {
        return const_reverse_iterator(begin());
    }

    const_reverse_iterator crend() const noexcept {
        return rend(); // Force const
    }

    reference front() noexcept {
        return *begin();
    }

    const_reference front() const noexcept {
        return *begin();
    }

    reference back() noexcept {
        return *rbegin();
    }

    const_reference back() const noexcept {
        return *rbegin();
    }

    Value* data() noexcept {
        return vals_.data();
    }

    const Value* data() const noexcept {
        return vals_.data();
    }

    template <std::integral Int>
    reference operator[](Int pos) noexcept {
        if constexpr (std::signed_integral<Int>)
            if (pos < 0)
                pos += size();
        return vals_[pos];
    }

    template <std::integral Int>
    const_reference operator[](Int pos) const noexcept {
        return const_cast<NdArray&>(*this)[pos];
    }

    template <std::integral Int>
    reference at(Int pos) {
        if constexpr (std::signed_integral<Int>)
            if (pos < 0)
                pos += size();
        if (pos < 0 or pos >= Int(size()))
            throw std::out_of_range(__func__);
        return vals_[pos];
    }

    template <std::integral Int>
    const_reference at(Int pos) const {
        return const_cast<NdArray&>(*this).at(pos);
    }

    void swap(NdArray& other) noexcept {
        dims_.swap(other.dims_);
        vals_.swap(other.vals_);
    }

    /** \} */

  public:
    /// \name Views
    /** \{ */

    template <size_t Rank>
    ArrayView<Value, Rank> view() {
        return {data(), sizes<Rank>()};
    }

    template <size_t Rank>
    ArrayView<const Value, Rank> view() const {
        return {data(), sizes<Rank>()};
    }

    template <size_t Rank>
    ArrayView<const Value, Rank> cview() const {
        return view();
    }

    template <size_t Rank>
    decltype(auto) operator[](ArrayIndex<Rank> k) {
        return view<Rank>()[k];
    }

    template <size_t Rank>
    decltype(auto) operator[](ArrayIndex<Rank> k) const {
        return view<Rank>()[k];
    }

    template <typename P, typename... Q>
    decltype(auto) operator()(P p, Q&&... q) {
        return view<1 + sizeof...(Q)>()(p, std::forward<Q>(q)...);
    }

    template <typename P, typename... Q>
    decltype(auto) operator()(P p, Q&&... q) const {
        return view<1 + sizeof...(Q)>()(p, std::forward<Q>(q)...);
    }

    /** \} */

  private:
    template <typename Other>
    using RebindVector = std::vector<
            Other,
            typename std::allocator_traits<Alloc>:: //
            template rebind_alloc<Other>>;

    RebindVector<ssize_t> dims_;
    RebindVector<Value> vals_;

  public:
    void serialize(auto& serializer) {
        serializer <=> dims_;
        serializer <=> vals_;
    }
};

} // namespace pre

namespace std {

template <typename... Args>
inline void swap(pre::NdArray<Args...>& lhs, pre::NdArray<Args...>& rhs) {
    lhs.swap(rhs);
}

} // namespace std
