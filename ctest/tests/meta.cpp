#include "../doctest.h"
#include <list>
#include <vector>
#include <pre/random>
#include <pre/meta>

TEST_CASE("Meta") {
    SUBCASE("IteratorRange") {
        pre::Pcg32 gen;
        std::list<int> vals0;
        std::vector<int> vals1;
        for (int j = 0; j < 256; j++) {
            vals0.emplace_back(gen());
            vals1.emplace_back(gen());
        }
        // Make ranges.
        auto range0 = pre::IteratorRange(vals0);
        auto range1 = pre::IteratorRange(vals1);
        auto itr0 = range0.begin();
        auto itr1 = range1.begin();
        CHECK(range0.size() == 256);
        CHECK(range1.size() == 256);
        // Also make ranges using iterators and pointers.
        auto also_range0 = pre::IteratorRange(vals0.begin(), vals0.end());
        auto also_range1 = pre::IteratorRange(&vals1[0], &vals1[0] + 256);
        auto also_itr0 = also_range0.begin();
        auto also_itr1 = also_range1.begin();
        CHECK(also_range0.size() == 256);
        CHECK(also_range1.size() == 256);
        for (int j = 0; j < 256; j++) {
            CHECK(*itr0 == *also_itr0);
            CHECK(*itr1 == *also_itr1);
            ++itr0;
            ++itr1;
            ++also_itr0;
            ++also_itr1;
        }
        CHECK(itr0 == range0.end());
        CHECK(itr1 == range1.end());
        CHECK(also_itr0 == also_range0.end());
        CHECK(also_itr1 == also_range1.end());
        // Indexing works even on list iterators, because the implementation
        // uses std::advance() to increment generically.
        int ind = gen(256);
        CHECK(range0[ind] == also_range0[ind]);
        CHECK(range1[ind] == also_range1[ind]);
    }
}
