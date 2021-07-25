OP1 = ['+', '-', '~', '!']
OP2 = ['+', '-', '*', '/', '%', '&', '|', '^', '>>', '<<']
#OPC = ['==', '!=', '<', '>', '<=', '>=', '&&', '||']

puts <<STR
namespace pre {

STR

for op1 in OP1
    puts <<STR
template <typename F, size_t R>
[[gnu::always_inline]]
constexpr auto operator#{op1}(const LazyArray<F, R>& arr) noexcept {
    return LazyArray([&](auto k) constexpr noexcept {
    return #{op1}arr(k); }, arr.sizes);
}

STR
end

for op2 in OP2
    puts <<STR
template <typename F, typename G, size_t R>
[[gnu::always_inline]]
constexpr auto operator#{op2}(
                    const LazyArray<F, R>& lhs, 
                    const LazyArray<G, R>& rhs) noexcept {
    ASSERT(lhs.sizes == rhs.sizes);
    return LazyArray([&](auto k) constexpr noexcept { 
    return lhs(k) #{op2} rhs(k); }, lhs.sizes);
}

STR
end

for op2 in OP2
    puts <<STR
template <typename F, concepts::arithmetic_or_complex T, size_t R>
[[gnu::always_inline]]
constexpr auto operator#{op2}(const LazyArray<F, R>& lhs, const T& rhs) noexcept {
    return LazyArray([&](auto k) constexpr noexcept { 
    return lhs(k) #{op2} rhs; }, lhs.sizes);
}

STR
end

for op2 in OP2
    puts <<STR
template <concepts::arithmetic_or_complex T, typename F, size_t R>
[[gnu::always_inline]]
constexpr auto operator#{op2}(const T& lhs, const LazyArray<F, R>& rhs) noexcept {
    return LazyArray([&](auto k) constexpr noexcept { 
    return lhs #{op2} rhs(k); }, rhs.sizes);
}

STR
end

puts <<STR

} // namespace pre

STR
