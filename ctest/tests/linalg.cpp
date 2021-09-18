#include "../doctest.h"
#include <pre/random>
#include <pre/Array>
#include <pre/math>

template <typename Field>
static bool is_approx_equal(Field lhs, Field rhs) {
    if constexpr (pre::concepts::matches<Field, std::complex>)
        return lhs.real() == Approx(rhs.real()).epsilon(1e-3) and
               lhs.imag() == Approx(rhs.imag()).epsilon(1e-3);
    else
        return lhs == Approx(rhs).epsilon(1e-3);
}

template <typename Field, size_t Rows, size_t... Cols>
static bool is_approx_equal(
    const pre::Array<Field, Rows, Cols...>& lhs,
    const pre::Array<Field, Rows, Cols...>& rhs) {
    for (size_t i = 0; i < Rows; i++)
        if (not is_approx_equal(lhs[i], rhs[i]))
            return false;
    return true;
}

template <typename Field, size_t Dim>
static bool is_approx_orthogonal(const pre::Array<Field, Dim, Dim>& x) {
    return is_approx_equal(
        pre::dot(x, pre::adjoint(x)), //
        pre::Array<Field, Dim, Dim>::identity());
}

template <typename Field, size_t Rows, size_t Cols>
pre::Array<Field, Rows, Cols> random_matrix(pre::Pcg32& gen) {
    pre::Array<Field, Rows, Cols> a;
    for (size_t i = 0; i < Rows; i++)
        for (size_t j = 0; j < Cols; j++)
            if constexpr (pre::concepts::matches<Field, std::complex>)
                a(i, j) = Field(gen(-1.0, 1.0), gen(-1.0, 1.0));
            else
                a(i, j) = Field(gen(-1.0, 1.0));
    return a;
}

template <typename Field, size_t Rows, size_t Cols>
void check_qr(pre::Pcg32& gen) {
    pre::Linalg<Field> linalg;
    pre::Array<Field, Rows, Cols> a = random_matrix<Field, Rows, Cols>(gen);
    pre::Array<Field, Rows, Cols> r = a;
    pre::Array<Field, Rows, Rows> q;
    CHECK_NOTHROW(linalg.qr(r, q));
    CHECK(is_approx_equal(pre::dot(q, r), a));
    CHECK(is_approx_orthogonal(q));
}

template <typename Field, size_t Dim>
void check_lu(pre::Pcg32& gen) {
    pre::Linalg<Field> linalg;
    // Do QR factorization to obtain orthogonal basis
    pre::Array<Field, Dim, Dim> a = random_matrix<Field, Dim, Dim>(gen);
    pre::Array<Field, Dim, Dim> q;
    CHECK_NOTHROW(linalg.qr(a, q));

    // Reset with random non-zero diagonal entries, expand with
    // randomized Q to obtain invertible matrix
    a = {};
    for (size_t k = 0; k < Dim; k++)
        if constexpr (pre::concepts::matches<Field, std::complex>)
            a(k, k) = Field(gen(0.1, 2.0), gen(0.1, 2.0));
        else
            a(k, k) = gen(0.1, 2.0);
    a = pre::dot(q, a, pre::adjoint(q));

    // Do factorization with and without pivoting
    pre::Array<Field, Dim, Dim> a0 = a;
    pre::Array<Field, Dim, Dim> a1 = a;
    pre::Array<int, Dim> p = {};
    CHECK_NOTHROW(linalg.lu(a0));
    CHECK_NOTHROW(linalg.lu(a1, p));

    // Solve linear system
    pre::Array<Field, Dim, 2> b0 = random_matrix<Field, Dim, 2>(gen), x0;
    pre::Array<Field, Dim, 2> b1 = random_matrix<Field, Dim, 2>(gen), x1;
    CHECK_NOTHROW(linalg.lu_solve(a0, b0, x0));
    CHECK_NOTHROW(linalg.lu_solve(a1, b1, x1, p));
    CHECK(is_approx_equal(pre::dot(a, x0), b0));
    CHECK(is_approx_equal(pre::dot(a, x1), b1));
}

template <typename Field, size_t Dim>
void check_chol(pre::Pcg32& gen) {
    pre::Linalg<Field> linalg;
    // Do QR factorization to obtain orthogonal basis
    pre::Array<Field, Dim, Dim> a = random_matrix<Field, Dim, Dim>(gen);
    pre::Array<Field, Dim, Dim> q;
    CHECK_NOTHROW(linalg.qr(a, q));

    // Reset with random positive diagonal entries, expand with
    // randomized Q to obtain hermitian positive definite matrix
    a = {};
    for (size_t k = 0; k < Dim; k++)
        a(k, k) = Field(gen(0.1, 2.0));
    a = pre::dot(q, a, pre::adjoint(q));
    for (size_t i = 0; i < Dim; i++) {
        a(i, i) = pre::real(a(i, i));
        for (size_t j = i + 1; j < Dim; j++) {
            a(j, i) = pre::conj(a(i, j));
        }
    }

    // Do factorization with and without pivoting
    pre::Array<Field, Dim, Dim> a0 = a;
    pre::Array<Field, Dim, Dim> a1 = a;
    pre::Array<int, Dim> p = {};
    CHECK_NOTHROW(linalg.chol(a0));
    CHECK_NOTHROW(linalg.chol(a1, p));
    auto p1 = pre::Array<Field, Dim, Dim>(pre::permutation_matrix(p));
    CHECK(is_approx_equal(pre::dot(pre::adjoint(a0), a0), a));
    CHECK(is_approx_equal(
        pre::dot(pre::adjoint(a1), a1), pre::dot(p1, a, pre::transpose(p1))));

    // Solve linear system
    pre::Array<Field, Dim, 2> b0 = random_matrix<Field, Dim, 2>(gen), x0;
    pre::Array<Field, Dim, 2> b1 = random_matrix<Field, Dim, 2>(gen), x1;
    CHECK_NOTHROW(linalg.chol_solve(a0, b0, x0));
    CHECK_NOTHROW(linalg.chol_solve(a1, b1, x1, p));
    CHECK(is_approx_equal(pre::dot(a, x0), b0));
    CHECK(is_approx_equal(pre::dot(a, x1), b1));
}

template <typename Field, size_t Rows, size_t Cols>
void check_svd(pre::Pcg32& gen) {
    pre::Linalg<Field> linalg;
    pre::Array<Field, Rows, Cols> a = random_matrix<Field, Rows, Cols>(gen);
    pre::Array<Field, Rows, Cols> x = a;
    pre::Array<Field, Rows, Rows> u;
    pre::Array<Field, Cols, Cols> v;
    CHECK_NOTHROW(linalg.svd(x, u, v));
    CHECK(is_approx_equal(pre::dot(u, x, v), a));
    CHECK(is_approx_orthogonal(u));
    CHECK(is_approx_orthogonal(v));
}

TEST_CASE_TEMPLATE("Linalg", Float, float, double) {
    pre::Pcg32 gen(getContextOptions()->rand_seed);
    SUBCASE("QR") {
        SUBCASE("Float") {
            check_qr<Float, 2, 2>(gen);
            check_qr<Float, 8, 3>(gen);
            check_qr<Float, 7, 8>(gen);
            check_qr<Float, 16, 22>(gen);
            check_qr<Float, 50, 50>(gen);
        }
        SUBCASE("std::complex<Float>") {
            check_qr<std::complex<Float>, 8, 3>(gen);
            check_qr<std::complex<Float>, 7, 8>(gen);
            check_qr<std::complex<Float>, 16, 22>(gen);
            check_qr<std::complex<Float>, 23, 11>(gen);
            check_qr<std::complex<Float>, 50, 50>(gen);
        }
    }
    SUBCASE("LU") {
        SUBCASE("Float") {
            check_lu<Float, 2>(gen);
            check_lu<Float, 9>(gen);
            check_lu<Float, 20>(gen);
            check_lu<Float, 50>(gen);
        }
        SUBCASE("std::complex<Float>") {
            check_lu<std::complex<Float>, 2>(gen);
            check_lu<std::complex<Float>, 9>(gen);
            check_lu<std::complex<Float>, 20>(gen);
            check_lu<std::complex<Float>, 50>(gen);
        }
    }
    SUBCASE("Cholesky") {
        SUBCASE("Float") {
            check_chol<Float, 4>(gen);
            check_chol<Float, 7>(gen);
            check_chol<Float, 20>(gen);
            check_chol<Float, 50>(gen);
        }
        SUBCASE("std::complex<Float>") {
            check_chol<std::complex<Float>, 4>(gen);
            check_chol<std::complex<Float>, 7>(gen);
            check_chol<std::complex<Float>, 20>(gen);
            check_chol<std::complex<Float>, 50>(gen);
        }
    }
    SUBCASE("SVD") {
        SUBCASE("Float") {
            check_svd<Float, 3, 3>(gen);
            check_svd<Float, 3, 5>(gen);
            check_svd<Float, 10, 11>(gen);
            check_svd<Float, 20, 13>(gen);
            check_svd<Float, 30, 35>(gen);
            check_svd<Float, 50, 50>(gen);
        }
        SUBCASE("std::complex<Float>") {
            check_svd<std::complex<Float>, 3, 3>(gen);
            check_svd<std::complex<Float>, 3, 5>(gen);
            check_svd<std::complex<Float>, 10, 11>(gen);
            check_svd<std::complex<Float>, 20, 13>(gen);
            check_svd<std::complex<Float>, 30, 35>(gen);
            check_svd<std::complex<Float>, 50, 50>(gen);
        }
    }
}
