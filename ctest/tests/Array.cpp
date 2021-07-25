#include "../doctest.h"
#include <sstream>
#include <pre/random>
#include <pre/Array>

TEST_CASE("Array") {
    pre::Pcg32 gen(getContextOptions()->rand_seed);

    SUBCASE("Rank-1") {
        auto arr = pre::Array{2, 3, 4, 7, 11};
        std::shuffle(arr.begin(), arr.end(), gen);
        SUBCASE("Rank and size") {
            CHECK(arr.rank() == 1);
            CHECK(arr.size() == 5);
            CHECK(arr.size() == arr.sizes()[0]);
        }
        SUBCASE("Utilities") {
            CHECK(arr.sum() == 2 + 3 + 4 + 7 + 11);
            CHECK(arr.prod() == 2 * 3 * 4 * 7 * 11);
            CHECK(arr[arr.argmin()] == 2);
            CHECK(arr[arr.argmax()] == 11);
            CHECK((arr.accum_sum().adjacent_difference() == arr).all());
        }
    }

    SUBCASE("Rank-3") {
        pre::Array<int, 2, 3, 4> arr;
        for (auto&& tmp0 : arr)
            for (auto&& tmp1 : tmp0)
                for (auto&& tmp2 : tmp1)
                    tmp2 = gen();
        CHECK(!(arr == 0).all());
        SUBCASE("Rank and size") {
            CHECK(arr.rank() == 3);
            CHECK(arr.sizes()[0] == 2);
            CHECK(arr.sizes()[1] == 3);
            CHECK(arr.sizes()[2] == 4);
            CHECK(arr.size() == arr.sizes()[0]);
            CHECK(arr[0].size() == arr.sizes()[1]);
            CHECK(arr[0][0].size() == arr.sizes()[2]);
        }
        SUBCASE("Flatten and reshape") {
            CHECK(arr.flatten().rank() == 1);
            CHECK((arr.flatten().reshape<2, 3, 4>() == arr).all());
        }
        SUBCASE("Stream operators") {
            pre::Array<int, 2, 3, 4> tmp;
            std::stringstream ss;
            ss << arr;
            ss >> tmp;
            CHECK((arr == tmp).all());
        }
    }

    SUBCASE("Vector/matrix") {
        SUBCASE("Dot vector with vector") {
            pre::Array<int, 3> lhs = {1, 2, 3};
            pre::Array<int, 3> rhs = {3, 5, 7};
            CHECK(pre::dot(lhs, rhs) == 1 * 3 + 2 * 5 + 3 * 7);
        }
        SUBCASE("Dot vector with matrix") {
            pre::Array<int, 2> lhs = {1, 2};
            pre::Array<int, 2, 4> rhs = {1, 2, 3, 4, 5, 6, 7, 8};
            CHECK((pre::dot(lhs, rhs) ==
                   pre::Array{
                           1 * 1 + 2 * 5, //
                           1 * 2 + 2 * 6, //
                           1 * 3 + 2 * 7, //
                           1 * 4 + 2 * 8})
                          .all());
        }
        SUBCASE("Dot matrix with vector") {
            pre::Array<int, 3, 2> lhs = {1, 2, 3, 4, 5, 6};
            pre::Array<int, 2> rhs = {1, 2};
            CHECK((pre::dot(lhs, rhs) ==
                   pre::Array{
                           1 * 1 + 2 * 2, //
                           3 * 1 + 4 * 2, //
                           5 * 1 + 6 * 2})
                          .all());
        }
        SUBCASE("Dot matrix with matrix") {
            pre::Array<int, 2, 2> lhs = {1, 2, 3, 4};
            pre::Array<int, 2, 3> rhs = {1, 2, 3, 4, 5, 6};
            CHECK((pre::dot(lhs, rhs) ==
                   pre::Array<int, 2, 3>{
                           1 * 1 + 2 * 4, 1 * 2 + 2 * 5, 1 * 3 + 2 * 6, //
                           3 * 1 + 4 * 4, 3 * 2 + 4 * 5, 3 * 3 + 4 * 6})
                          .all());
        }
    }
}

TEST_CASE("ArrayView") {
    pre::Pcg32 gen(getContextOptions()->rand_seed);

    SUBCASE("Rank and size") {
        pre::Array<int, 5, 8, 7, 2> arr;
        CHECK(arr->rank() == 4);
        CHECK(arr->size() == arr->sizes[0]);
        CHECK(arr->sizes[0] == 5);
        CHECK(arr->sizes[1] == 8);
        CHECK(arr->sizes[2] == 7);
        CHECK(arr->sizes[3] == 2);
    }

    SUBCASE("Rank-1 slicing") {
        pre::Array<int, 7> arr = {1, 2, 3, 4, 5, 6, 7};
        SUBCASE("Slice") {
            auto tmp = arr(pre::Slice(1, 4));
            CHECK(tmp.size() == 3);
            CHECK(tmp[0] == 2);
            CHECK(tmp[1] == 3);
            CHECK(tmp[2] == 4);
        }
        SUBCASE("Slice reversal") {
            auto tmp = arr(pre::Slice(4, 1));
            CHECK(tmp.size() == 3);
            CHECK(tmp[0] == 4);
            CHECK(tmp[1] == 3);
            CHECK(tmp[2] == 2);
        }
        SUBCASE("Slice to end") {
            auto tmp = arr(pre::Slice(2, -1));
            CHECK(tmp.size() == 5);
            CHECK(tmp[0] == 3);
            CHECK(tmp[1] == 4);
            CHECK(tmp[2] == 5);
            CHECK(tmp[3] == 6);
            CHECK(tmp[4] == 7);
        }
    }

    SUBCASE("Rank-2 slicing") {
        pre::Array<int, 12, 18> arr;
        for (auto&& tmp0 : arr)
            for (auto&& tmp1 : tmp0)
                tmp1 = gen();
        CHECK(!(arr == 0).all());

        SUBCASE("Rank-1 slice") {
            auto tmp = arr(pre::Slice(2, 9), 8);
            CHECK(tmp.rank() == 1);
            CHECK(tmp.size() == 7);
            for (int i = 0; i < 7; i++)
                CHECK(&tmp[i] == &arr[i + 2][8]);
        }
        SUBCASE("Rank-2 slice") {
            auto tmp = arr(pre::Slice(2, 9), pre::Slice(7, 18));
            CHECK(tmp.rank() == 2);
            CHECK(tmp.sizes[0] == 7);
            CHECK(tmp.sizes[1] == 11);
            for (int i = 0; i < 7; i++)
                for (int j = 0; j < 11; j++)
                    CHECK(&tmp(i, j) == &arr[i + 2][j + 7]);
        }
        SUBCASE("Empty rank-1 slice") {
            auto tmp = arr(pre::Slice(2, 2), 11);
            CHECK(tmp.rank() == 1);
            CHECK(tmp.size() == 0);
        }
        SUBCASE("Empty rank-2 slice") {
            auto tmp = arr(pre::Slice(2, 2), pre::Slice(0, -1));
            CHECK(tmp.rank() == 2);
            CHECK(tmp.sizes[0] == 0);
            CHECK(tmp.sizes[1] == 18);
            CHECK(tmp.empty());
        }
        SUBCASE("Matrix-only methods") {
            CHECK(arr->rows() == arr->sizes[0]);
            CHECK(arr->cols() == arr->sizes[1]);
            CHECK(arr->row(7).size() == 18);
            CHECK(arr->col(3).size() == 12);

            // Diagonal.
            auto tmp = arr->diag();
            CHECK(tmp.size() == 12);
            for (int k = 0; k < 12; k++)
                CHECK(&tmp[k] == &arr[k][k]);

            // Above diagonal.
            tmp = arr->diag(2);
            CHECK(tmp.rank() == 1);
            CHECK(tmp.size() == 12);
            for (int k = 0; k < 12; k++)
                CHECK(&tmp[k] == &arr[k][k + 2]);

            // Transpose.
            CHECK(arr->transpose().rows() == 18);
            CHECK(arr->transpose().cols() == 12);
            CHECK(&arr->transpose()(8, 4) == &arr[4][8]);
        }
    }
}
