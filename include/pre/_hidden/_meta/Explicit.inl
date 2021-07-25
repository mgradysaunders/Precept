/*-*- C++ -*-*/
#pragma once

namespace pre {

template <typename T>
struct Explicit {
  public:
    template <typename... Args>
    constexpr void construct(Args&&... args) {
        new (&bytes[0]) T(std::forward<Args>(args)...);
    }

    template <typename... Args>
    constexpr void reconstruct(Args&&... args) {
        deconstruct();
        construct(std::forward<Args>(args)...);
    }

    constexpr void deconstruct() noexcept {
        if constexpr (!concepts::trivially_destructible<T>)
            (**this).~T();
    }

  public:
    constexpr T& operator*() noexcept {
        return *reinterpret_cast<T*>(&bytes[0]);
    }

    constexpr const T& operator*() const noexcept {
        return *reinterpret_cast<const T*>(&bytes[0]);
    }

    constexpr T* operator->() noexcept {
        return reinterpret_cast<T*>(&bytes[0]);
    }

    constexpr const T* operator->() const noexcept {
        return reinterpret_cast<T*>(&bytes[0]);
    }

  public:
    alignas(alignof(T)) std::byte bytes[sizeof(T)];
};

} // namespace pre
