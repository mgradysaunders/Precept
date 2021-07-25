#include "../doctest.h"
#include <pre/random>
#include <pre/Array>
#include <pre/math>

template <typename Float, size_t N>
static bool is_approx(
        const pre::Array<Float, N>& lhs, //
        const pre::Array<Float, N>& rhs) {
    for (size_t k = 0; k < N; k++)
        if (!(lhs[k] == Approx(rhs[k]).epsilon(0.05)))
            return false;
    return true;
}

TEST_CASE("Math") {
    pre::Pcg32 gen(getContextOptions()->rand_seed);
    SUBCASE("min, max, clamp, nthpow") {
        CHECK(pre::min(2, 3) == 2);
        CHECK(pre::max(2, 3) == 3);
        CHECK(pre::clamp(0.5, 0.6, 0.8) == 0.6);
        CHECK(pre::clamp(0.7, 0.6, 0.8) == 0.7);
        CHECK(pre::clamp(0.9, 0.6, 0.8) == 0.8);
        CHECK(pre::nthpow(2.0, 4) == 16.0);
        CHECK(pre::nthpow(4.0, -2) == 1.0 / 16.0);
    }
    SUBCASE("repeat, mirror, ispow2, roundpow2") {
        CHECK(pre::repeat(3, 7) == 3);
        CHECK(pre::repeat(10, 7) == 3);
        CHECK(pre::repeat(-1, 7) == 6);
        CHECK(pre::mirror(3, 7) == 3);
        CHECK(pre::mirror(11, 7) == 2);
        CHECK(pre::mirror(-1, 7) == 0);
        CHECK(pre::ispow2(1 << 17));
        CHECK(pre::ispow2(1 << 17) - 1 == false);
        CHECK(pre::roundpow2(0) == 1);
        CHECK(pre::roundpow2(7) == 8);
        CHECK(pre::roundpow2(500) == 512);
    }
    SUBCASE("first1, rotl, rotr, bit_swap, bit_reverse") {
        CHECK(pre::first1(1) == 0);
        CHECK(pre::first1(1 << 10) == 10);
        CHECK(std::uint32_t(0x22334411UL) ==
              pre::rotl(std::uint32_t(0x11223344UL), 8));
        CHECK(std::uint32_t(0x44112233UL) ==
              pre::rotr(std::uint32_t(0x11223344UL), 8));
        CHECK(0b0010UL == pre::bit_swap(0b1000UL, 1, 3));
        std::uint32_t bits0 = 0b11101010'10001111'10110111'10001010UL;
        std::uint32_t bits1 = 0b01010001'11101101'11110001'01010111UL;
        CHECK(bits0 == pre::bit_reverse(bits1));
    }
    SUBCASE("gcd_bezout, gcd") {
        for (int k = 0; k < 1024; k++) {
            int a = gen() % (1 << 12);
            int b = gen() % (1 << 12);
            int x = 0, y = 0;
            int d = pre::gcd_bezout(a, b, &x, &y);
            CHECK(d == pre::gcd(a, b));
            CHECK(d == a * x + b * y); // Bezout condition
        }
    }
}

TEST_CASE_TEMPLATE("Math", Float, float, double) {
    pre::Pcg32 gen(getContextOptions()->rand_seed);
    SUBCASE("float_incr, float_decr") {
        for (int k = 0; k < 1024; k++) {
            Float x = pre::generate_canonical<Float>(gen) * 128 - 64;
            CHECK(pre::float_incr(x) ==
                  std::nextafter(x, +std::numeric_limits<Float>::infinity()));
            CHECK(pre::float_decr(x) ==
                  std::nextafter(x, -std::numeric_limits<Float>::infinity()));
        }
        // No increment past infinity.
        CHECK(pre::float_incr(+pre::numeric_limits<Float>::infinity()) ==
              +pre::numeric_limits<Float>::infinity());
        CHECK(pre::float_decr(-pre::numeric_limits<Float>::infinity()) ==
              -pre::numeric_limits<Float>::infinity());
    }
    SUBCASE("fast_floor, fast_ceil, fast_round, fast_trunc") {
        for (int k = 0; k < 1024; k++) {
            Float x = pre::generate_canonical<Float>(gen) * 128 - 64;
            if ((x - pre::floor(x)) == Float(0.5))
                continue;
            CHECK(Float(pre::fast_floor(x)) == pre::floor(x));
            CHECK(Float(pre::fast_ceil(x)) == pre::ceil(x));
            CHECK(Float(pre::fast_round(x)) == pre::round(x));
            CHECK(Float(pre::fast_trunc(x)) == pre::trunc(x));
        }
    }
    SUBCASE("sinpi, cospi, sinpi_min, cospi_min") {
        for (int k = 0; k < 1024; k++) {
            Float x = pre::generate_canonical<Float>(gen) * 128 - 64;
            Float y1 = pre::sinpi(x);
            Float z1 = pre::cospi(x);
            Float y2 = pre::sin(pre::numeric_constants<Float>::M_pi() * x);
            Float z2 = pre::cos(pre::numeric_constants<Float>::M_pi() * x);
            CHECK(y1 == Approx(y2).epsilon(0.01));
            CHECK(z1 == Approx(z2).epsilon(0.01));
        }
        for (int j = 0; j < 256; j++) {
            Float x1 = pre::generate_canonical<Float>(gen) * 128 - 64;
            Float x2 = pre::generate_canonical<Float>(gen) * 2 + x1;
            Float ymin = pre::sinpi_min(x1, x2) - Float(1e-4);
            Float ymax = pre::sinpi_max(x1, x2) + Float(1e-4);
            for (int k = 0; k < 1024; k++) {
                Float t = Float(k) / 1024;
                Float x = pre::lerp(t, x1, x2);
                Float y = pre::sinpi(x);
                CHECK(ymin < y);
                CHECK(ymax > y);
            }
        }
    }
    SUBCASE("erfinv") {
        for (int j = 0; j < 4096; j++) {
            Float x = pre::generate_canonical<Float>(gen) * 32 - 16;
            Float y = pre::erf(x);
            if (pre::abs(y) < Float(0.999))
                CHECK(pre::erfinv(y) == Approx(x).epsilon(0.05));
        }
    }
    pre::Mat3<Float> rgb_to_xyz = pre::rgb_to_xyz(
            pre::Array{Float(0.7350), Float(0.2650)},
            pre::Array{Float(0.2740), Float(0.7170)},
            pre::Array{Float(0.1670), Float(0.0090)},
            pre::Array{Float(1), Float(1), Float(1)});
    for (int check = 0; check < 512; check++) {
        pre::Vec3<Float> rgb = {
                pre::generate_canonical<Float>(gen),
                pre::generate_canonical<Float>(gen),
                pre::generate_canonical<Float>(gen)};
        pre::Vec3<Float> xyz = pre::rgb_to_xyz(rgb);
        SUBCASE("srgbenc, srgbdec") {
            CHECK(is_approx(rgb, pre::srgbdec(pre::srgbenc(rgb))));
        }
        SUBCASE("rgb_to_xyz, xyz_to_rgb, xyz_to_xyy, xyy_to_xyz") {
            CHECK(is_approx(rgb, pre::xyz_to_rgb(xyz)));
            CHECK(is_approx(xyz, pre::dot(rgb_to_xyz, rgb)));
            CHECK(is_approx(xyz, pre::xyy_to_xyz(pre::xyz_to_xyy(xyz))));
            CHECK(is_approx(xyz, pre::lab_to_xyz(pre::xyz_to_lab(xyz))));
        }
        SUBCASE("xyz_to_luv, luv_to_xyz") {
            pre::Vec3<Float> w = {
                    pre::generate_canonical<Float>(gen),
                    pre::generate_canonical<Float>(gen),
                    pre::generate_canonical<Float>(gen)};
            CHECK(is_approx(xyz, pre::luv_to_xyz(w, pre::xyz_to_luv(w, xyz))));
        }
    }
}
