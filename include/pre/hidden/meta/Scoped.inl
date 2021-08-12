/*-*- C++ -*-*/
#pragma once

namespace pre {

/// An RAII scoped management helper.
///
/// \tparam F  A function object to call on construction.
/// \tparam G  A function object to call on destruction.
///
template <typename F, typename G>
struct Scoped {
    constexpr Scoped(F&& f, G&& g)
        : f_(std::forward<F>(f)), g_(std::forward<G>(g)) {
        std::invoke(f_);
    }
    constexpr ~Scoped() {
        std::invoke(g_);
    }

    // Non-copyable.
    Scoped(const Scoped&) = delete;

  private:
    F f_;
    G g_;
};

template <typename F, typename G>
Scoped(F&&, G&&) -> Scoped<F, G>;

} // namespace pre
