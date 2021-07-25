/*-*- C++ -*-*/
#pragma once

namespace pre {

namespace concepts {

/// Is type a pointer?
template <typename T>
concept pointer = std::is_pointer_v<T>;

/// Is type a reference?
template <typename T>
concept reference = std::is_reference_v<T>;

/// Is type an lvalue reference?
template <typename T>
concept lvalue_reference = std::is_lvalue_reference_v<T>;

/// Is type an rvalue reference?
template <typename T>
concept rvalue_reference = std::is_rvalue_reference_v<T>;

/// Is type a subclass of X?
template <typename T, typename X>
concept subclass = std::is_base_of_v<X, T>;

/// Is type a superclass of X?
template <typename T, typename X>
concept superclass = std::is_base_of_v<T, X>;

// Is instantiation?
template <typename, template <typename...> typename>
struct is_instantiation : std::false_type {};

// Is instantiation? true specialization.
template <typename... Ts, template <typename...> typename X>
struct is_instantiation<X<Ts...>, X> : std::true_type {};

/// Matches template?
template <typename T, template <typename...> typename X>
concept matches = is_instantiation<T, X>::value;

/// Does not match template?
template <typename T, template <typename...> typename X>
concept not_matches = !is_instantiation<T, X>::value;

/// Matches `std::basic_istream`?
template <typename T>
concept istream = requires {
    requires subclass<
            T, std::basic_istream<
                       typename T::char_type, typename T::traits_type>>;
};

/// Matches `std::basic_ostream`?
template <typename T>
concept ostream = requires {
    requires subclass<
            T, std::basic_ostream<
                       typename T::char_type, typename T::traits_type>>;
};

// Is arithmetic?
template <typename T>
concept arithmetic = std::is_arithmetic_v<T>;

// Is arithmetic or enum?
template <typename T>
concept arithmetic_or_enum = std::is_arithmetic_v<T> || std::is_enum_v<T>;

/// Is arithmetic or standard complex?
template <typename T>
concept arithmetic_or_complex = arithmetic<T> || matches<T, std::complex>;

/// Is floating point or standard complex?
template <typename T>
concept floating_point_or_complex =
        std::floating_point<T> || matches<T, std::complex>;

// Is functor?
template <typename T>
concept functor = requires {
    requires std::is_member_function_pointer_v<decltype(&T::operator())>;
};

/// Is minimal container?
template <typename T>
concept minimal_container = requires(T& c) {
    requires std::default_initializable<T>;
    requires std::destructible<T>;
    requires std::signed_integral<typename T::difference_type>;
    requires std::unsigned_integral<typename T::size_type>;
    // clang-format off
    { c.begin()  } -> std::convertible_to<typename T::iterator>;
    { c.end()    } -> std::convertible_to<typename T::iterator>;
    { c.cbegin() } -> std::convertible_to<typename T::const_iterator>;
    { c.cend()   } -> std::convertible_to<typename T::const_iterator>;
    { c.size()   } -> std::convertible_to<typename T::size_type>;
    // clang-format on
};

/// Is allocator aware?
template <typename T>
concept allocator_aware = requires(T& c) {
    requires std::constructible_from<T, const T&, typename T::allocator_type>;
    requires std::constructible_from<T, T&&, typename T::allocator_type>;
    requires std::assignable_from<T, const T&>;
    requires std::assignable_from<T, T&&>;
    { c.get_allocator() }
    ->std::convertible_to<typename T::allocator_type>;
};

/// Is sequence constructible?
template <typename T>
concept sequence_constructible = requires {
    requires std::constructible_from<
            T, const typename T::value_type*, const typename T::value_type*>;
};

/// Is trivially_destructible?
template <typename T>
concept trivially_destructible = std::is_trivially_destructible_v<T>;

/// Is random access range of given value type?
template <typename T, typename Value>
concept random_access_range_of = requires {
    requires std::ranges::random_access_range<T>;
    requires std::convertible_to<std::ranges::range_value_t<T>*, Value*>;
};

} // namespace concepts

} // namespace pre
