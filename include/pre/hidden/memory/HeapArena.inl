/*-*- C++ -*-*/
#pragma once

namespace pre {

/// A heap-allocated memory arena.
template <typename Alloc = std::allocator<std::byte>>
class HeapArena {
  public:
    typedef Alloc allocator_type;

    typedef std::allocator_traits<Alloc> allocator_traits;

    HeapArena(size_t block_size = 0, const Alloc& alloc = Alloc()) noexcept
        : block_size_(block_size), //
          free_blocks_(alloc),     //
          full_blocks_(alloc), alloc_(alloc) {
        // Round up to 256 byte interval.
        block_size_ += 255u;
        block_size_ &= ~255u;
        if (block_size_ == 0)
            block_size_ = 65536;
        // Allocate initial block, reserve blocks.
        block_.size = block_size_;
        block_.begin = alloc_.allocate(block_.size);
        block_.offset = 0;
        free_blocks_.reserve(4);
        full_blocks_.reserve(4);
    }

    HeapArena(const HeapArena&) = delete;

    HeapArena(HeapArena&& other) noexcept
        : block_size_(steal(other.block_size_)), //
          block_(steal(other.block_)),
          free_blocks_(std::move(other.free_blocks_)),
          full_blocks_(std::move(other.full_blocks_)),
          alloc_(std::move(other.alloc_)) {
    }

    HeapArena(HeapArena&& other, const Alloc& alloc) : alloc_(alloc) {
        // We need to be able to deallocate the pointers we steal!
        if (alloc_ != other.alloc)
            throw std::invalid_argument(__func__);
        block_size_ = steal(other.block_size_);
        block_ = steal(other.block_);
        free_blocks_ = std::move(other.free_blocks_);
        full_blocks_ = std::move(other.full_blocks_);
    }

    ~HeapArena() {
        alloc_.deallocate(block_.begin, block_.size);
        for (Block& block : free_blocks_)
            alloc_.deallocate(block.begin, block.size);
        for (Block& block : full_blocks_)
            alloc_.deallocate(block.begin, block.size);
    }

    HeapArena& operator=(const HeapArena&) = delete;

    HeapArena& operator=(HeapArena&& other) noexcept {
        block_size_ = steal(other.block_size_);
        block_ = steal(other.block_);
        free_blocks_ = std::move(other.free_blocks_);
        full_blocks_ = std::move(other.full_blocks_);
        if constexpr (allocator_traits::
                              propagate_on_container_move_assignment::value)
            alloc_ = std::move(other.alloc_);
        return *this;
    }

  public:
    /// Allocate bytes.
    void* allocate(size_t sz) {
        // Round up to 16 byte interval.
        sz += 15u;
        sz &= ~15u;
        if (sz == 0)
            return nullptr;

        if (block_.size < block_.offset + sz) {
            full_blocks_.emplace_back(block_);
            if (free_blocks_.size() == 0 || free_blocks_.back().size < sz) {
                // Allocate block.
                block_.size = std::max(block_size_, sz);
                block_.begin = alloc_.allocate(block_.size);
                block_.offset = 0;
            }
            else {
                // Use free block.
                block_ = free_blocks_.back();
                free_blocks_.pop_back();
            }
        }
        std::byte* pos = block_.begin + block_.offset;
        block_.offset += sz;
        return static_cast<void*>(pos);
    }

    /// Allocate given type.
    template <typename T>
    T* allocate(size_t count = 1) {
        return static_cast<T*>(allocate(sizeof(T) * count));
    }

    /// Clear.
    void clear() {
        // Clear current block.
        block_.offset = 0;
        // Convert full blocks to free blocks.
        free_blocks_.reserve(free_blocks_.size() + full_blocks_.size());
        for (Block& block : full_blocks_) {
            free_blocks_.push_back(block);
            free_blocks_.back().offset = 0;
        }
        full_blocks_.clear();
    }

    /// Clear and deallocate.
    void reset() {
        block_.offset = 0;
        for (Block& block : free_blocks_)
            alloc_.deallocate(block.begin, block.size);
        for (Block& block : full_blocks_)
            alloc_.deallocate(block.begin, block.size);
        free_blocks_.clear();
        full_blocks_.clear();
    }

    void swap(HeapArena& other) {
        if (this != &other) {
            std::swap(block_size_, other.block_size_);
            std::swap(block_, other.block_);
            std::swap(free_blocks_, other.free_blocks_);
            std::swap(full_blocks_, other.full_blocks_);
            if constexpr (allocator_traits::propagate_on_container_swap::value)
                std::swap(alloc_, other.alloc_);
        }
    }

  private:
    /// A memory block.
    struct Block {
        std::byte* begin; ///< Pointer to bytes.
        size_t offset;    ///< Offset.
        size_t size;      ///< Size.
    };

    size_t block_size_;

    Block block_;

    template <typename T>
    using RebindAlloc = typename allocator_traits::template rebind_alloc<T>;

    std::vector<Block, RebindAlloc<Block>> free_blocks_;

    std::vector<Block, RebindAlloc<Block>> full_blocks_;

    RebindAlloc<std::byte> alloc_;
};

} // namespace pre

template <typename Alloc>
inline void* operator new(size_t sz, pre::HeapArena<Alloc>& arena) {
    return arena.allocate(sz);
}

template <typename Alloc>
inline void* operator new[](size_t sz, pre::HeapArena<Alloc>& arena) {
    return arena.allocate(sz);
}

template <typename Alloc>
inline void operator delete(void*, pre::HeapArena<Alloc>&) noexcept {
    // Do nothing
}

template <typename Alloc>
inline void operator delete[](void*, pre::HeapArena<Alloc>&) noexcept {
    // Do nothing
}
