/*-*- C++ -*-*/
#pragma once

namespace pre {

/// A static stack.
///
/// A static stack (first-in last-out), where elements are
/// pushed onto and popped from the top/back.
///
template <typename Value, size_t MaxSize>
struct StaticStack : ArrayLike<StaticStack<Value, MaxSize>, Value> {
  public:
    constexpr StaticStack() noexcept = default;

    template <std::input_or_output_iterator Iterator>
    constexpr StaticStack(Iterator from, Iterator to) {
        while (from != to)
            push(*from++);
    }

  public:

    /// \name Container API
    /** \{ */

    constexpr size_t size() const noexcept {
        return top_;
    }

    constexpr size_t max_size() const noexcept {
        return MaxSize;
    }

    constexpr size_t capacity() const noexcept {
        return MaxSize;
    }

    constexpr Value* begin() noexcept {
        return &values_[0];
    }

    constexpr const Value* begin() const noexcept {
        return &values_[0];
    }

    constexpr Value* end() noexcept {
        return &values_[0] + top_;
    }

    constexpr const Value* end() const noexcept {
        return &values_[0] + top_;
    }

    constexpr bool empty() const noexcept {
        return top_ == 0;
    }

    constexpr bool full() const noexcept {
        return top_ == MaxSize;
    }

    constexpr void clear() noexcept {
        for (auto& value : *this)
            value.~Value();
        top_ = 0;
    }

    /** \} */

public:
    /// \name Stack API
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

    /// Pop and return top/back value.
    ///
    /// \throw std::runtime_error  If empty.
    ///
    constexpr Value pop() {
        if (top_ == 0)
            throw std::runtime_error(__func__);
        Value& value = values_[--top_];
        Value result = std::move(value);
        value.~Value();
        return result;
    }

    /// Push top/back value.
    ///
    /// \throw std::length_error  If full.
    ///
    constexpr void push(const Value& value) {
        if (top_ >= MaxSize)
            throw std::length_error(__func__);
        values_[top_++] = value;
    }

    /** \} */

  public:
    void serialize(auto& serializer) {
        serializer <=> values_;
        serializer <=> top_;
    }

  private:
    Value values_[MaxSize] = {};

    size_t top_ = 0;
};

} // namespace pre
