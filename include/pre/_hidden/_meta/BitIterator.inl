/*-*- C++ -*-*/
#pragma once

namespace pre {

template <std::integral Int>
class BitIterator {
  public:
    using value_type = int;
    using reference = const int&;
    using difference_type = int;
    constexpr BitIterator() noexcept = default;
    constexpr BitIterator(Int bits, int bit = 0) noexcept
        : bits_(bits), bit_(bit) {
        while (bit_ < max_bit and (bits_ & (Int(1) << bit_)) == 0)
            bit_++;
    }
    constexpr const int& operator*() const noexcept {
        return bit_;
    }
    constexpr BitIterator& operator++() noexcept {
        while (bit_ < max_bit) {
            bit_++;
            if (bit_ == max_bit or (bits_ & (Int(1) << bit_)))
                break;
        }
        return *this;
    }
    constexpr BitIterator operator++(int) noexcept {
        BitIterator copy = *this;
        operator++();
        return copy;
    }
    constexpr auto operator<=>(const BitIterator& other) const noexcept {
        return bit_ <=> other.bit_;
    }
    constexpr bool operator==(const BitIterator& other) const noexcept {
        return bit_ == other.bit_;
    }
    constexpr bool operator!=(const BitIterator& other) const noexcept {
        return bit_ != other.bit_;
    }

  private:
    static constexpr int max_bit = sizeof(Int) * 8;
    Int bits_ = ~Int(0);
    int bit_ = 0;
};

} // namespace pre
