/*-*- C++ -*-*/
#pragma once

namespace pre {

/// A reference counter for a reference countable pointer.
template <typename T>
class RefPtr {
  public:
    constexpr RefPtr() noexcept = default;

    constexpr RefPtr(std::nullptr_t) noexcept {
    }

    RefPtr(const RefPtr& other) noexcept : ptr_(other.ptr_) {
        incr_ref(ptr_);
    }

    RefPtr(RefPtr&& other) noexcept : ptr_(other.ptr_) {
        other.ptr_ = nullptr;
    }

    explicit RefPtr(T* ptr) : ptr_(ptr) {
        incr_ref(ptr_);
    }

    /// Steal pointer, dummy int for overload
    explicit RefPtr(T* ptr, int) : ptr_(ptr) {
    }

    ~RefPtr() {
        decr_ref(ptr_);
    }

    RefPtr& operator=(const RefPtr& other) noexcept {
        if (ptr_ != other.ptr_) {
            decr_ref(ptr_);
            ptr_ = other.ptr_;
            incr_ref(ptr_);
        }
        return *this;
    }

    RefPtr& operator=(RefPtr&& other) noexcept {
        decr_ref(ptr_);
        ptr_ = other.ptr_;
        other.ptr_ = nullptr;
        return *this;
    }

  public:
    T* get() noexcept {
        return ptr_;
    }

    T* get() const noexcept {
        return ptr_;
    }

    void reset() noexcept {
        decr_ref(ptr_);
        ptr_ = nullptr;
    }

    template <typename U>
    void reset(U* ptr) noexcept {
        reset();
        ptr_ = ptr;
        incr_ref(ptr_);
    }

    void swap(RefPtr& other) noexcept {
        std::swap(ptr_, other.ptr_);
    }

    long use_count() const noexcept {
        return ptr_ ? long(ptr_->ref_count) : 0;
    }

  public:
    auto operator<=>(const RefPtr&) const = default;

    bool operator==(const RefPtr&) const = default;

    decltype(auto) operator*() noexcept {
        return (*ptr_);
    }

    decltype(auto) operator*() const noexcept {
        return (*ptr_);
    }

    decltype(auto) operator->() noexcept {
        return ptr_;
    }

    decltype(auto) operator->() const noexcept {
        return ptr_;
    }

    operator bool() const noexcept {
        return ptr_ != nullptr;
    }

    template <typename U>
    operator RefPtr<U>&() noexcept requires std::is_base_of_v<U, T> {
        return reinterpret_cast<RefPtr<U>&>(*this);
    }

    template <typename U>
    operator const RefPtr<U>&() const noexcept
            requires std::is_base_of_v<U, T> {
        return reinterpret_cast<const RefPtr<U>&>(*this);
    }

  private:
    T* ptr_ = nullptr;
};

template <concepts::subclass<RefCountable> T>
RefPtr(T*) -> RefPtr<T>;

template <typename U, typename T>
inline auto ref_cast(const RefPtr<T>& ref) {
    return RefPtr(reinterpret_cast<U*>(ref.get()));
}

template <typename T, typename... Args>
inline auto make_ref(Args&&... args) {
    return RefPtr(new T(std::forward<Args>(args)...));
}

} // namespace pre
