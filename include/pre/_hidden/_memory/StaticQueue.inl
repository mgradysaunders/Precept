/*-*- C++ -*-*/
#pragma once

namespace pre {

/// A static queue.
///
/// A static queue (first-in first-out), where elements are
/// pushed onto the top/back and popped from the bottom/front.
///
/// \note
/// The implementation minimizes move operations by tracking
/// the effective bottom and top of the queue. That is, popping
/// shifts the bottom up instead of moving the elements
/// down. When pushing shifts the top past the queue capacity, the
/// implementation moves all elements down and resets the bottom
/// to zero.
///
template <typename Value, size_t MaxSize>
struct StaticQueue : ArrayLike<StaticQueue<Value, MaxSize>, Value> {
  public:
    // Sanity check.
    static_assert(MaxSize > 0);

  public:
    constexpr StaticQueue() = default;

    template <std::input_or_output_iterator Iterator>
    constexpr StaticQueue(Iterator from, Iterator to) {
        while (from != to)
            push(*from++);
    }

  public:
    /// \name Container API
    /** \{ */

    constexpr size_t size() const noexcept {
        return top_ - bot_;
    }

    constexpr size_t max_size() const noexcept {
        return MaxSize;
    }

    constexpr size_t capacity() const noexcept {
        return MaxSize;
    }

    constexpr void clear() noexcept {
        for (auto& value : *this)
            value.~Value();
        bot_ = 0;
        top_ = 0;
    }

    constexpr bool empty() const noexcept {
        return bot_ == top_;
    }

    constexpr bool full() const noexcept {
        return bot_ == 0 && top_ == MaxSize;
    }

    constexpr Value* begin() noexcept {
        return &values_[0] + bot_;
    }

    constexpr const Value* begin() const noexcept {
        return &values_[0] + bot_;
    }

    constexpr Value* end() noexcept {
        return &values_[0] + top_;
    }

    constexpr const Value* end() const noexcept {
        return &values_[0] + top_;
    }

    /** \} */

public:

    /// \name Queue
    /** \{ */

    /// Bottom/front value with empty check.
    ///
    /// \throw std::runtime_error  If empty.
    ///
    constexpr decltype(auto) bottom() {
        if (empty())
            throw std::runtime_error(__func__);
        return this->front();
    }

    /// Bottom/front value with empty check, const variant.
    ///
    /// \throw std::runtime_error  If empty.
    ///
    constexpr decltype(auto) bottom() const {
        if (empty())
            throw std::runtime_error(__func__);
        return this->front();
    }

    /// Top/back value with empty check.
    ///
    /// \throw std::runtime_error  If empty.
    ///
    constexpr decltype(auto) top() {
        if (empty())
            throw std::runtime_error(__func__);
        return this->back();
    }

    /// Top/back value with empty check, const variant.
    ///
    /// \throw std::runtime_error  If empty.
    ///
    constexpr decltype(auto) top() const {
        if (empty())
            throw std::runtime_error(__func__);
        return this->back();
    }

    /// Pop and return bottom/front value.
    ///
    /// \throw std::runtime_error  If empty.
    ///
    constexpr Value pop() {
        if (bot_ == top_)
            throw std::runtime_error(__func__);
        Value& value = values_[bot_++];
        Value res = std::move(value);
        value.~Value();
        if (bot_ == top_) // Empty now?
            bot_ = top_ = 0;
        return res;
    }

    /// Push top/back value.
    ///
    /// \throw std::length_error  If full.
    ///
    constexpr void push(const Value& value) {
        if (top_ >= MaxSize) {
            if (bot_ == 0)
                throw std::length_error(__func__);
            Value* itr0 = &values_[0];
            Value* itr1 = &values_[0] + bot_;
            for (; itr1 < &values_[0] + top_; ++itr0, ++itr1) {
                *itr0 = std::move(*itr1);
                itr1->~Value();
            }
            top_ -= bot_;
            bot_ = 0;
        }
        values_[top_++] = value;
    }

    /** \} */

  public:
    void serialize(auto& serializer) {
        serializer <=> values_;
        serializer <=> bot_;
        serializer <=> top_;
    }

  private:
    Value values_[MaxSize] = {};

    size_t bot_ = 0;

    size_t top_ = 0;
};

} // namespace pre
