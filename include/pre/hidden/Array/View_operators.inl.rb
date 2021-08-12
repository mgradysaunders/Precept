OP2 = ['+', '-', '*', '/', '%', '&', '|', '^', '>>', '<<']

puts <<STR
namespace pre {

STR

for op2 in OP2
    puts <<STR
template <typename T, typename U, size_t R>
[[gnu::always_inline]]
constexpr ArrayView<T, R> operator#{op2}=(
          ArrayView<T, R> lhs, ArrayView<U, R> rhs) noexcept
{
    ASSERT(lhs.sizes == rhs.sizes);
    auto itrlhs = lhs.begin();
    auto itrrhs = rhs.begin();
    for (; itrrhs != rhs.end(); ++itrlhs, ++itrrhs)
        *itrlhs #{op2}= *itrrhs;
    return lhs;
}

STR
end

for op2 in OP2
    puts <<STR
template <typename T, size_t R>
[[gnu::always_inline]]
constexpr ArrayView<T, R> operator#{op2}=(
          ArrayView<T, R> lhs, const std::type_identity_t<T>& rhs) noexcept
{
    for (auto itr = lhs.begin(); itr != lhs.end(); ++itr)
        *itr #{op2}= rhs;
    return lhs;
}

STR
end

puts <<STR

} // namespace pre

STR
