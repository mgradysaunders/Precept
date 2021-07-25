#include "../doctest.h"
#include <pre/memory>

TEST_CASE("StaticStack") {
    pre::StaticStack<int, 8> stack;

    for (int k = 0; k < 8; k++) {
        // Should not throw yet.
        CHECK_NOTHROW(stack.push(k));
        // Top == most recently pushed.
        CHECK(stack.top() == k);
        // Top == back.
        CHECK(stack.top() == stack.back());
    }

    // Front == bottom.
    CHECK(stack.front() == stack.bottom());
    // Should now be full.
    CHECK(stack.full());
    // If full, push should throw.
    CHECK_THROWS(stack.push(8));

    for (int k = 0; k < 8; k++) {
        // Element at k should be k.
        CHECK(stack[k] == k);
        // Reverse indexing should also work.
        CHECK(stack[-(k + 1)] == 7 - k);
    }

    for (int k = 7; k >= 0; k--) {
        // Top == most recently pushed.
        CHECK(stack.top() == k);
        // Pop should be most recently pushed.
        CHECK(stack.pop() == k);
    }

    // Should now be empty.
    CHECK(stack.empty());
    // If empty, bottom should throw.
    CHECK_THROWS(stack.bottom());
    // If empty, top should throw.
    CHECK_THROWS(stack.top());
    // If empty, front should still not throw.
    CHECK_NOTHROW(stack.front());
    // If empty, back should still not throw.
    CHECK_NOTHROW(stack.back());
    // If empty, pop should throw.
    CHECK_THROWS(stack.pop());
}
