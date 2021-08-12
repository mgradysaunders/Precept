/*-*- C++ -*-*/
#pragma once

namespace pre {

/// A heap-allocated memory stack.
///
/// This structure is meant to optimize heap-allocations in contexts
/// where memory in use can be checkpointed and rolled back using stack-style
/// push and pop operations.
///
template <typename Alloc = std::allocator<std::byte>>
struct HeapStack {
  public:
    typedef Alloc allocator_type;

    typedef std::allocator_traits<Alloc> allocator_traits;

    HeapStack(size_t block_size = 0, const Alloc& alloc = Alloc())
        : block_size_(block_size), alloc_(alloc), pushes_(alloc) {
        if (block_size_ == 0)
            block_size_ = 65536;
        block_tail_ = block_allocate_(0);
        pushes_.reserve(8);
    }

    HeapStack(const HeapStack&) = delete;

    HeapStack(HeapStack&& other) noexcept
        : block_size_(other.block_size_),
          block_tail_(steal(other.block_tail_)),
          alloc_(std::move(other.alloc_)), pushes_(std::move(other.pushes_)) {
    }

    ~HeapStack() {
        clear();
        block_deallocate_(block_tail_);
        block_tail_ = nullptr;
    }

    HeapStack& operator=(const HeapStack&) = delete;

    HeapStack& operator=(HeapStack&& other) noexcept {
        this->~HeapStack();
        block_size_ = other.block_size_;
        block_tail_ = steal(other.block_tail_);
        alloc_ = std::move(other.alloc_);
        pushes_ = std::move(other.pushes_);
        return *this;
    }

  public:
    void* allocate(size_t n) {
        n = (n + 15U) & ~15U;
        ensure_tail_can_allocate_(n);
        void* ptr = block_tail_->top;
        block_tail_->top += n;
        return ptr;
    }

    void clear() noexcept {
        // Find head.
        Block* head = block_tail_;
        while (head->prev)
            head = head->prev;

        // Deallocate all remaining blocks.
        Block* block = head->next;
        while (block) {
            Block* next = block->next;
            block_deallocate_(block);
            block = next;
        }

        // Clear.
        block_tail_ = head;
        block_tail_->next = nullptr;
        block_tail_->top = block_tail_->begin;
        pushes_.clear();
    }

    /// Create scoped push object to automatically pop on destruction.
    auto scoped_push() {
        return Scoped([&]() { push(); }, [&]() { pop(); });
    }

    /// Push current state.
    void push() {
        pushes_.push_back({block_tail_, block_tail_->top});
    }

    /// Pop state.
    void pop() {
        // Tail block.
        block_tail_ = pushes_.back().first;
        block_tail_->top = pushes_.back().second;
        pushes_.pop_back();

        // Clear all next blocks.
        Block* block = block_tail_->next;
        while (block) {
            block->top = block->begin;
            block = block->next;
        }
    }

  private:
    /// A block.
    struct alignas(16) Block {
        Block* next;      ///< Pointer to next block.
        Block* prev;      ///< Pointer to previous block.
        std::byte* begin; ///< Pointer to stack memory begin.
        std::byte* end;   ///< Pointer to stack memory end.
        std::byte* top;   ///< Pointer to stack memory top.
    };

    size_t block_size_ = 65536;

    Block* block_tail_ = nullptr;

    template <typename T>
    using RebindAlloc =
            typename std::allocator_traits<Alloc>::template rebind_alloc<T>;

    template <typename T>
    using RebindVector = std::vector<T, RebindAlloc<T>>;

    RebindAlloc<std::byte> alloc_;

    RebindVector<std::pair<Block*, std::byte*>> pushes_;

  private:
    Block* block_allocate_(size_t least_n = 0) {
        if (least_n < block_size_)
            least_n = block_size_;
        auto* ptr = alloc_.allocate(sizeof(Block) + least_n);
        auto* block = reinterpret_cast<Block*>(ptr);
        block->next = nullptr;
        block->prev = nullptr;
        block->begin = ptr + sizeof(Block);
        block->end = block->begin + least_n;
        block->top = block->begin;
        return block;
    }

    void block_deallocate_(Block* block) {
        alloc_.deallocate(
                reinterpret_cast<std::byte*>(block),
                sizeof(Block) + block->end - block->begin);
    }

    void ensure_tail_can_allocate_(size_t n) {
        if (block_tail_->top + n > block_tail_->end) {
            if (block_tail_->next) {
                block_tail_ = block_tail_->next;
                ensure_tail_can_allocate_(n); // Recurse.
            }
            else {
                block_tail_->next = block_allocate_(n);
                block_tail_->next->prev = block_tail_;
                block_tail_ = block_tail_->next;
            }
        }
    }
};

} // namespace pre
