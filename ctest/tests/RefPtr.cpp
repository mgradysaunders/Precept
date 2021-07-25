#include "../doctest.h"
#include <pre/random>
#include <pre/memory>

TEST_CASE("RefPtr") {
    // Call counters
    static int Foo_ctor_calls = 0;
    static int Foo_dtor_calls = 0;

    class Foo final : virtual public pre::RefCountable {
      public:
        Foo(int v = 0) : value(v) {
            ++Foo_ctor_calls;
        }
        ~Foo() {
            ++Foo_dtor_calls;
        }

      public:
        int value = 0;
    };

    pre::Pcg32 gen(getContextOptions()->rand_seed);
    {
        // Allocate and move around
        std::vector<pre::RefPtr<Foo>> foos;
        for (int repeat = 0; repeat < 8; repeat++) {
            for (int iter = 0; iter < 64; iter++) {
                foos.emplace_back(new Foo(gen()));
            }
            for (int iter = 0; iter < 8; iter++) {
                auto pos = gen(foos.size());
                auto foo = foos[pos];
                CHECK(foo.use_count() == 2);
                foos.erase(foos.begin() + pos);
                CHECK(foo.use_count() == 1);
                foos.push_back(foo);
            }
            std::sort(
                    foos.begin(), foos.end(),
                    [](const auto& lhs, const auto& rhs) {
                        return lhs->value < rhs->value;
                    });
        }
    }
    CHECK(Foo_ctor_calls == Foo_dtor_calls);
}
