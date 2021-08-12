puts <<STR
namespace pre {

STR

puts <<STR
template <typename T>
using VecView = ArrayView<T, 1>;

STR

for m in [1, 2, 3, 4]
    puts <<STR
template <typename T>
using Vec#{m} = Array<T, #{m}>;

STR
end

puts <<STR
template <typename T>
using MatView = ArrayView<T, 2>;

STR

for m in [1, 2, 3, 4]
    puts <<STR
template <typename T>
using Mat#{m} = Array<T, #{m}, #{m}>;

STR
end

for m in [1, 2, 3, 4]
for n in [1, 2, 3, 4]
    puts <<STR
template <typename T>
using Mat#{m}x#{n} = Array<T, #{m}, #{n}>;

STR
end
end

puts <<STR

} // namespace pre

STR
