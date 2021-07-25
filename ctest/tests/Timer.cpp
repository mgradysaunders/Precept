#include "../doctest.h"
#include <thread>
#include <pre/Timer>

TEST_CASE("Timer") {
    pre::SystemTimer timer;
    SUBCASE("Read after sleep") {
        std::this_thread::sleep_for(std::chrono::seconds(1));
        CHECK(timer.nanoseconds() > std::int64_t(0.999 * 1e9));
        CHECK(timer.microseconds() > std::int64_t(0.999 * 1e6));
        CHECK(timer.milliseconds() > std::int64_t(0.999 * 1e3));
        CHECK(timer.seconds() > 0.999);
    }
    SUBCASE("Read again after reset") {
        timer.reset();
        CHECK(timer.nanoseconds() < std::int64_t(0.999 * 1e9));
        CHECK(timer.microseconds() < std::int64_t(0.999 * 1e6));
        CHECK(timer.milliseconds() < std::int64_t(0.999 * 1e3));
        CHECK(timer.seconds() < 0.999);
    }
}
