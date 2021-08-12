/*-*- C++ -*-*/
#pragma once

namespace pre {

/// An atomic signed lock-free integer.
///
/// \note
/// This is provided by the C++20 standard as `std::atomic_signed_lock_free`,
/// however it is not quite available in GCC yet. This declaration should be
/// removed once it becomes available.
///
using atomic_signed_lock_free = std::conditional_t<
        std::atomic_ptrdiff_t::is_always_lock_free,
        std::atomic_ptrdiff_t,
        std::conditional_t<
                std::atomic_long::is_always_lock_free,
                std::atomic_long,
                std::atomic_int>>;

/// A reference countable interface.
class RefCountable {
  public:
    virtual ~RefCountable() = default;

  public:
    atomic_signed_lock_free ref_count = 0;
};

inline void incr_ref(RefCountable* ptr) noexcept {
    if (ptr)
        ++ptr->ref_count;
}

inline void decr_ref(RefCountable* ptr) noexcept {
    if (ptr) {
        if (--ptr->ref_count == 0) {
            delete ptr;
        }
    }
}

} // namespace pre
