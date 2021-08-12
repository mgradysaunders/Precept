OP1 = ['+', '-', '~', '!']
OP2 = ['+', '-', '*', '/', '%', '&', '|', '^', '>>', '<<']
OPC = ['==', '!=', '<', '>', '<=', '>=', '&&', '||']

puts <<STR
namespace pre {

STR

#------------------------------------------------------------------------------

for op1 in OP1
    puts <<STR
template <typename T, size_t... N>
[[gnu::always_inline]]
constexpr auto operator#{op1}(const Array<T, N...>& arr) noexcept
{
    using U = decltype(#{op1}T());
    Array<U, N...> res;
    auto itrarr = arr.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrarr, ++itrres)
        *itrres = #{op1}*itrarr;
    return res;
}

STR
end

#------------------------------------------------------------------------------

for op2 in (OP2 + OPC)
    puts <<STR
template <typename T, typename U, size_t... N>
[[gnu::always_inline]]
constexpr auto operator#{op2}(
                    const Array<T, N...>& lhs, 
                    const Array<U, N...>& rhs) noexcept
{
    using V = decltype(T() #{op2} U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrrhs, ++itrres)
        *itrres = *itrlhs #{op2} *itrrhs;
    return res;
}

STR
end

for op2 in OP2
    puts <<STR
template <typename T, typename U, size_t... N>
[[gnu::always_inline]]
constexpr Array<T, N...>& operator#{op2}=(
          Array<T, N...>& lhs, const Array<U, N...>& rhs) noexcept
{
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs #{op2}= *itrrhs;
    return lhs;
}

STR
end

#------------------------------------------------------------------------------

for op2 in (OP2 + OPC)
    puts <<STR
template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]]
constexpr auto operator#{op2}(const Array<T, N...>& lhs, const U& rhs) noexcept
{
    using V = decltype(T() #{op2} U());
    Array<V, N...> res;
    auto itrlhs = lhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrlhs, ++itrres)
        *itrres = *itrlhs #{op2} rhs;
    return res;
}

STR
end

for op2 in OP2
    puts <<STR
template <typename T, concepts::not_array U, size_t... N>
[[gnu::always_inline]]
constexpr Array<T, N...>& operator#{op2}=(Array<T, N...>& lhs, const U& rhs) noexcept
{
    for (auto itrlhs = lhs.begin(); itrlhs != lhs.end(); ++itrlhs)
        *itrlhs #{op2}= rhs;
    return lhs;
}

STR
end

#------------------------------------------------------------------------------

for op2 in (OP2 + OPC)
    puts <<STR
template <concepts::not_array T, typename U, size_t... N>
[[gnu::always_inline]]
constexpr auto operator#{op2}(const T& lhs, const Array<U, N...>& rhs) noexcept
{
    using V = decltype(T() #{op2} U());
    Array<V, N...> res;
    auto itrrhs = rhs.begin();
    auto itrres = res.begin();
    for (; itrres != res.end(); ++itrrhs, ++itrres)
        *itrres = lhs #{op2} *itrrhs;
    return res;
}

STR
end

puts <<STR

} // namespace pre

STR
