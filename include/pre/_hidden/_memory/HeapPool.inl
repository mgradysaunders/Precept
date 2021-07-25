/*-*- C++ -*-*/
#pragma once

namespace pre {

/// A heap-allocated memory pool.
template <typename Alloc = std::allocator<std::byte>>
class HeapPool {
  public:
    typedef Alloc allocator_type;

    typedef std::allocator_traits<Alloc> allocator_traits;

    /// Constructor.
    ///
    /// \param[in] elem_size
    /// Element size in bytes.
    ///
    /// \param[in] pool_size
    /// Elements per pool. If left at default value of 0, this
    /// is chosen so that the size of each pool is approximately 64KB.
    ///
    /// \param[in] alloc
    /// Allocator.
    ///
    HeapPool(
            size_t elem_size,
            size_t pool_size = 0,
            const Alloc& alloc = Alloc()) noexcept
        : elem_size_(elem_size), //
          pool_size_(pool_size), //
          alloc_(alloc) {
        if (elem_size_ < sizeof(void*))
            elem_size_ = sizeof(void*);
        if (pool_size_ == 0) {
            pool_size_ = 65536 / elem_size;
            if (pool_size_ == 0)
                pool_size_ = 1;
        }
    }

    HeapPool(const HeapPool&) = delete;

    HeapPool(HeapPool&& other) noexcept
        : elem_size_(steal(other.elem_size_)),
          pool_size_(steal(other.pool_size_)), //
          first_(steal(other.first_)), alloc_(std::move(other.alloc_)) {
    }

    HeapPool(HeapPool&& other, const Alloc& alloc) : alloc_(alloc) {
        // We need to be able to deallocate the pointers we steal!
        if (alloc_ == other.alloc)
            throw std::invalid_argument(__func__);
        elem_size_ = steal(other.elem_size_);
        pool_size_ = steal(other.pool_size_);
        first_ = steal(other.first_);
    }

    ~HeapPool() {
        reset();
    }

    HeapPool& operator=(const HeapPool&) = delete;

    HeapPool& operator=(HeapPool&& other) noexcept {
        elem_size_ = steal(other.elem_size_);
        pool_size_ = steal(other.pool_size_);
        first_ = steal(other.first_);
        if constexpr (allocator_traits::
                              propagate_on_container_move_assignment::value)
            alloc_ = std::move(other.alloc_);
        return *this;
    }

  public:
    [[nodiscard]] void* allocate() {
        Pool* pool = first_;
        while (pool && pool->first_free == nullptr)
            pool = pool->next;
        if (pool == nullptr)
            pool = pool_allocate_();
        std::byte* elem = pool->first_free;
        std::memcpy(&pool->first_free, pool->first_free, sizeof(void*));
        pool->count++;
        return elem;
    }

    void deallocate(void* ptr) {
        if (ptr == nullptr)
            return;
        std::byte* elem = static_cast<std::byte*>(ptr);
        Pool* pool = first_;
        while (pool) {
            // Element in range of pool?
            std::ptrdiff_t pool_size = elem_size_ * pool_size_;
            std::ptrdiff_t elem_diff = elem - pool->begin;
            if (elem_diff >= 0 && elem_diff < pool_size) {
                // Element in range of pool, but pointer is not
                // aligned to element boundary?
                if (elem_diff % std::ptrdiff_t(elem_size_) != 0)
                    throw std::logic_error(__func__); // Garbage pointer!
                break;
            }
            else {
                pool = pool->next;
            }
        }
        // Element not in range of any pool?
        if (!pool)
            throw std::logic_error(__func__); // Garbage pointer!
        // Prepend.
        std::memcpy(elem, &pool->first_free, sizeof(void*));
        pool->first_free = elem;
        pool->count--;
    }

    /// Clear.
    void clear() noexcept {
        for (Pool* pool = first_; pool; pool = pool->next)
            pool_clear_(pool);
    }

    /// Clear and deallocate.
    void reset() noexcept {
        for (Pool* pool = first_; pool;) {
            Pool* next = pool->next;
            pool_deallocate_(pool);
            pool = next;
        }
        first_ = nullptr;
    }

    void swap(HeapPool& other) noexcept {
        if (this != &other) {
            std::swap(elem_size_, other.elem_size_);
            std::swap(pool_size_, other.pool_size_);
            std::swap(first_, other.first_);
            if constexpr (allocator_traits::propagate_on_container_swap::value)
                std::swap(alloc_, other.alloc_);
        }
    }

  public:
    /// \name Observers
    /** \{ */

    size_t bytes_per_elem() const noexcept {
        return elem_size_;
    }

    size_t elems_per_pool() const noexcept {
        return pool_size_;
    }

    Alloc get_allocator() {
        return alloc_;
    }

    /** \} */

  protected:
    /// Pool type.
    struct alignas(16) Pool {
        Pool* next;            ///< Pointer to next pool.
        std::byte* begin;      ///< Pointer to elements.
        std::byte* first_free; ///< Pointer to first free element.
        size_t count;          ///< Element allocation count.
    };

    size_t elem_size_ = 0; ///< Element size in bytes.
    size_t pool_size_ = 0; ///< Pool size in elements.

    /// First pool in list.
    Pool* first_ = nullptr;

    /// Allocator.
    typename allocator_traits::template rebind_alloc<std::byte> alloc_;

  private:
    Pool* pool_allocate_() {
        auto* ptr = alloc_.allocate(sizeof(Pool) + elem_size_ * pool_size_);
        Pool* pool = reinterpret_cast<Pool*>(ptr);
        pool->next = nullptr;
        pool->begin = ptr + sizeof(Pool);
        pool->count = 0;
        pool_clear_(pool);
        // Append.
        if (first_ == nullptr) {
            first_ = pool; // Set first.
        }
        else {
            Pool* last = first_;
            while (last->next)
                last = last->next;
            last->next = pool;
        }
        return pool;
    }

    void pool_deallocate_(Pool* pool) {
        alloc_.deallocate(
                reinterpret_cast<std::byte*>(pool),
                sizeof(Pool) + elem_size_ * pool_size_);
    }

    void pool_clear_(Pool* pool) {
        // Link all elements sequentially.
        pool->first_free = pool->begin;
        for (size_t index = 0; index + 1 < pool_size_; index++) {
            std::byte* elem0 = pool->begin + (index + 0) * elem_size_;
            std::byte* elem1 = pool->begin + (index + 1) * elem_size_;
            std::memcpy(elem0, &elem1, sizeof(void*));
        }
        std::byte* elem0 = pool->begin + (pool_size_ - 1) * elem_size_;
        std::byte* elem1 = nullptr;
        std::memcpy(elem0, &elem1, sizeof(void*));
        // Zero count.
        pool->count = 0;
    }
};

/// A heap-allocated pool for a given object.
template <typename Obj, typename Alloc = std::allocator<std::byte>>
class ObjectHeapPool final : public HeapPool<Alloc> {
  public:
    using Base = HeapPool<Alloc>;

    ObjectHeapPool(size_t pool_size, const Alloc& alloc = Alloc())
        : Base(sizeof(Obj), pool_size, alloc) {
    }

    ObjectHeapPool(ObjectHeapPool&& other) noexcept : Base(other) {
    }

    ObjectHeapPool(ObjectHeapPool&& other, const Alloc& alloc)
        : Base(other, alloc) {
    }

    ~ObjectHeapPool() = default;

    Obj* allocate() {
        return static_cast<Obj*>(static_cast<Base&>(*this).allocate());
    }

    template <typename... Args>
    Obj* create(Args&&... args) {
        return new (allocate()) Obj(std::forward<Args>(args)...);
    }

    void destroy(Obj* object) {
        if (object) {
            object->~Obj();
            static_cast<Base&>(*this).deallocate(object);
        }
    }
};

} // namespace pre
