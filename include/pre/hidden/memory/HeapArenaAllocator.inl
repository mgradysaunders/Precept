/*-*- C++ -*-*/
#pragma once

namespace pre {

/// A standard-compatible memory arena allocator.
template <typename T, typename Alloc = std::allocator<std::byte>>
class HeapArenaAllocator {
  public:
    typedef T value_type;

    typedef std::true_type propagate_on_container_copy_assignment;

    typedef std::true_type propagate_on_container_move_assignment;

    typedef std::true_type propagate_on_container_swap;

    typedef std::false_type is_always_equal;

  public:
    HeapArenaAllocator(size_t block_size = 0, const Alloc& alloc = Alloc())
        : arena_(new HeapArena<Alloc>(block_size, alloc)) {
    }

    template <typename U>
    HeapArenaAllocator(const HeapArenaAllocator<U, Alloc>& other)
        : arena_(other.arena_) {
    }

    template <typename U>
    HeapArenaAllocator(HeapArenaAllocator<U, Alloc>&& other)
        : arena_(std::move(other.arena_)) {
    }

    template <typename U>
    HeapArenaAllocator& operator=(const HeapArenaAllocator<U, Alloc>& other) {
        if (this != &other) {
            this->arena_ = other.arena_;
        }
        return *this;
    }

    template <typename U>
    HeapArenaAllocator& operator=(HeapArenaAllocator<U, Alloc>&& other) {
        this->arena_ = std::move(other.arena_);
        return *this;
    }

    void clear() {
        arena_->clear();
    }

    void reset() {
        arena_->reset();
    }

    [[nodiscard]] T* allocate(size_t n) {
        return static_cast<T*>(arena_->allocate(sizeof(T) * n));
    }

    void deallocate(T*, size_t) {
    }

    template <typename U>
    bool operator==(const HeapArenaAllocator<U, Alloc>& other) const {
        return arena_.get() == other.arena_.get();
    }

    template <typename U>
    bool operator!=(const HeapArenaAllocator<U, Alloc>& other) const {
        return arena_.get() != other.arena_.get();
    }

  private:
    std::shared_ptr<HeapArena<Alloc>> arena_;

    template <typename, typename>
    friend class HeapArenaAllocator;
};

} // namespace pre
