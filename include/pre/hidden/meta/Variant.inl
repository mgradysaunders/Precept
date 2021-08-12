/*-*- C++ -*-*/
#pragma once

namespace pre {

template <typename... Types>
class Variant {
  public:
    constexpr Variant() = default;
    constexpr Variant(const Variant&) = default;
    constexpr Variant(Variant&&) = default;
    template <typename Type>
    constexpr Variant(Type&& v) : var_(std::forward<Type>(v)) {
    }
    template <typename... Other>
    constexpr explicit Variant(Other&&... args)
        : var_(std::forward<Other>(args)...) {
    }

  public:
    template <typename Type>
    constexpr bool is() const noexcept {
        return std::holds_alternative<Type>(var_);
    }
    template <typename Type>
    constexpr decltype(auto) get() {
        return std::get<Type>(var_);
    }
    template <typename Type>
    constexpr decltype(auto) get() const {
        return std::get<Type>(var_);
    }
    template <size_t Index>
    constexpr decltype(auto) get() {
        return std::get<Index>(var_);
    }
    template <size_t Index>
    constexpr decltype(auto) get() const {
        return std::get<Index>(var_);
    }
    constexpr decltype(auto) index() const noexcept {
        return var_.index();
    }
    constexpr decltype(auto) valueless_by_exception() const noexcept {
        return var_.valueless_by_exception();
    }
    template <typename Functor>
    constexpr decltype(auto) dispatch(Functor&& func) {
        return std::visit(std::forward<Functor>(func), var_);
    }
    template <typename Functor>
    constexpr decltype(auto) dispatch(Functor&& func) const {
        return std::visit(std::forward<Functor>(func), var_);
    }

  public:
    constexpr Variant& operator=(const Variant&) = default;
    constexpr Variant& operator=(Variant&&) = default;
    template <typename Type>
    constexpr Variant& operator=(Type&& v) {
        var_ = std::forward<Type>(v);
        return *this;
    }
    constexpr auto operator<=>(const Variant& other) const {
        return var_ <=> other.var_;
    }
    constexpr bool operator==(const Variant& other) const {
        return var_ == other.var_;
    }
    constexpr bool operator!=(const Variant& other) const {
        return var_ != other.var_;
    }
    constexpr operator std::variant<Types...>&() noexcept {
        return var_;
    }
    constexpr operator const std::variant<Types...>&() const noexcept {
        return var_;
    }
  public:
    void serialize(auto& serializer) {
        serializer <=> var_;
    }

  private:
    std::variant<Types...> var_;
};

} // namespace pre
