#include "../doctest.h"
#include <pre/memory>

TEST_CASE("StaticQueue") {
    pre::StaticQueue<int, 8> queue;

    for (int k = 0; k < 8; k++) {
        // Should not throw yet.
        CHECK_NOTHROW(queue.push(k));
        // Top == most recently pushed.
        CHECK(queue.top() == k);
        // Top == back.
        CHECK(queue.top() == queue.back());
    }

    // Should now be full.
    CHECK(queue.full());
    // If full, push should throw.
    CHECK_THROWS(queue.push(8));

    for (int k = 0; k < 8; k++) {
        // Element at k should be k.
        CHECK(queue[k] == k);
        // Reverse indexing should also work.
        CHECK(queue[-(k + 1)] == 7 - k);
    }

    for (int k = 0; k < 8; k++) {
        // Front == bottom.
        CHECK(queue.front() == queue.bottom());
        // Front == least recently pushed.
        CHECK(queue.front() == k);
        // Pop should be least recently pushed.
        CHECK(queue.pop() == k);
    }

    // Should now be empty.
    CHECK(queue.empty());
    // If empty, bottom should throw.
    CHECK_THROWS(queue.bottom());
    // If empty, top should throw.
    CHECK_THROWS(queue.top());
    // If empty, front should still not throw.
    CHECK_NOTHROW(queue.front());
    // If empty, back should still not throw.
    CHECK_NOTHROW(queue.back());
    // If empty, pop should throw.
    CHECK_THROWS(queue.pop());
}
