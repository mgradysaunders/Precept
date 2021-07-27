/*-*- C++ -*-*/
#pragma once

namespace pre {

template <
        concepts::floating_point_or_complex Field,
        typename Alloc = std::allocator<Field>>
struct Linalg {
  public:
    using Float = to_floating_point_t<Field>;
    using FloatLimits = numeric_limits<Float>;

  private:
    /// A cache to prevent excessive heap allocations.
    struct Cache {
        auto scoped_push() {
            return stack_.scoped_push();
        }
        void clear() noexcept {
            stack_.clear();
        }
        template <typename T>
        VecView<T> vector(ssize_t n) {
            T* ptr = static_cast<T*>(stack_.allocate(sizeof(T) * n));
            return {ptr, n};
        }
        template <typename T>
        MatView<T> matrix(ssize_t m, ssize_t n) {
            T* ptr = static_cast<T*>(stack_.allocate(sizeof(T) * m * n));
            const ssize_t sizes[] = {m, n};
            return {ptr, &sizes[0]};
        }

      private:
        HeapStack<Alloc> stack_;
    };

    Cache cache_;

  private:
    static void expect(bool cond, const char* what) {
        if (!cond)
            throw std::runtime_error(what);
    }

  public:
    /// \name Miscellaneous
    /** \{ */

    static constexpr void load_identity(MatView<Field> a) noexcept {
        a = Field(0);
        int m = a.rows();
        int n = a.cols();
        for (int k = 0; k < pre::min(m, n); k++)
            a(k, k) = Field(1);
    }

    static constexpr void load_iota(VecView<int> v) noexcept {
        int n = 0;
        for (int& i : v)
            i = n++;
    }

    static constexpr void force_upper_triangular(MatView<Field> a) noexcept {
        int m = a.rows();
        int n = a.cols();
        for (int i = 1; i < m; i++)
            a(i, Slice(0, pre::min(i, n))) = Field(0);
    }

    static constexpr void force_lower_triangular(MatView<Field> a) noexcept {
        int m = a.rows();
        int n = a.cols();
        for (int j = 1; j < n; j++)
            a(Slice(0, pre::min(j, m)), j) = Field(0);
    }

    static constexpr void do_adjoint(MatView<Field> a) noexcept {
        ASSERT(a.is_square());
        int m = a.rows();
        for (int i = 0; i < m; i++) {
            a(i, i) = pre::conj(a(i, i));
            for (int j = i + 1; j < m; j++) {
                auto a0 = a(i, j);
                auto a1 = a(j, i);
                a(i, j) = pre::conj(a1);
                a(j, i) = pre::conj(a0);
            }
        }
    }

    static constexpr int find_pivot(VecView<Field> v) noexcept {
        Float tmp_max = 0;
        int piv = 0;
        int n = v.size();
        for (int k = 0; k < n; k++) {
            Float tmp = pre::abs(v[k]);
            if (tmp_max < tmp) {
                tmp_max = tmp;
                piv = k;
            }
        }
        return piv;
    }

    /** \} */

  public:
    /// \name Norms
    /** \{ */

    Float norm1(VecView<const Field> v) {
        return pre::sum(pre::abs(*v));
    }

    Float norm2(VecView<const Field> v) {
        int n = v.size();
        switch (n) {
        case 0: return 0;
        case 1: return std::abs(v[0]);
        case 2: return std::hypot(std::abs(v[0]), std::abs(v[1]));
        default: break;
        }
        auto _ = cache_.scoped_push();
        auto absv = cache_.template vector<Float>(n);
        absv = pre::abs(*v);
        auto absv_max = *std::max_element(absv.begin(), absv.end());
        if (absv_max == 0)
            return 0;
        else if (
                absv_max * absv_max * n > FloatLimits::max() ||
                absv_max < FloatLimits::min_squarable()) {
            absv /= absv_max;
            return pre::sqrt(pre::sum(pre::norm(*absv))) * absv_max;
        }
        else {
            return pre::sqrt(pre::sum(pre::norm(*absv)));
        }
    }

    Float normalize1(VecView<Field> v) {
        if (auto len = norm1(v); len > 0) {
            v /= len;
            return len;
        }
        else {
            return 0;
        }
    }

    Float normalize2(VecView<Field> v) {
        if (auto len = norm2(v); len > 0) {
            v /= len;
            return len;
        }
        else {
            return 0;
        }
    }

    /** \} */

  public:
    /// \name Factorization: LU
    /** \{ */

    /// Factorize.
    ///
    /// \param[inout] a  Matrix to factorize.
    /// \param[out]   p  _Optional_. Permutation vector.
    ///
    /// \returns
    /// Returns sign of pivot permutation, either `+1` or `-1`. If no
    /// pivoting, this is always `+1`.
    ///
    static int lu(MatView<Field> a, VecView<int> p = {}) {
        ASSERT(a.is_square());
        ASSERT(!p || a.size() == p.size());
        int n = a.size(); // Dimension.
        int s = 1;        // Permutation sign.
        if (p)
            load_iota(p);
        for (int j = 0; j < n; j++) {
            if (p) { // Optionally pivot.
                int piv = find_pivot(a(Slice(j, n), j)) + j;
                if (piv != j) {
                    std::swap(p[j], p[piv]);
                    a.swap_rows(j, piv);
                    s = -s;
                }
            }
            expect(pre::abs(a(j, j)) > FloatLimits::min_invertible(),
                   "non-singular matrix");
            // Factorize.
            Field fac = Float(1) / a(j, j);
            for (int i = j + 1; i < n; i++) {
                a(i, j) *= fac;
                for (int k = j + 1; k < n; k++)
                    a(i, k) -= a(i, j) * a(j, k);
            }
        }
        return s;
    }

    /// Use factorization to solve linear system \f$ AX = B \f$.
    ///
    /// \param[in]  a  Matrix of coefficients, factorized in-place.
    /// \param[in]  b  Matrix of coefficients.
    /// \param[out] x  Matrix of solutions.
    /// \param[in]  p  _Optional_. Permutation vector.
    ///
    static void lu_solve(
            MatView<const Field> a,
            MatView<const Field> b,
            MatView<Field> x,
            VecView<const int> p = {}) {
        ASSERT(a.is_square());
        ASSERT(b.rows() == a.rows() && //
               x.rows() == a.cols() && //
               x.cols() == b.cols());
        ASSERT(!p || a.size() == p.size());
        int m = b.rows();
        int n = b.cols();
        for (int j = 0; j < n; j++) {
            // Solve Ly = b.
            for (int i = 0; i < m; i++) {
                auto k = Slice(0, i);
                x(i, j) = b(p ? p[i] : i, j) - dot(*a(i, k), *x(k, j));
            }
            // Solve Ux = y.
            for (int i = m - 1; i >= 0; i--) {
                auto k = Slice(i + 1, m);
                x(i, j) -= dot(*a(i, k), *x(k, j));
                x(i, j) /= a(i, i);
            }
        }
    }

    /** \} */

  public:
    /// \name Factorization: Cholesky
    /** \{ */

    /// Factorize.
    ///
    /// \param[inout] a  Matrix t to factorize.
    /// \param[out]   p  _Optional_. Permutation vector. 
    ///
    static void chol(MatView<Field> a, VecView<int> p = {}) {
        ASSERT(a.is_square());
        ASSERT(!p || a.size() == p.size());
        int n = a.size();
        if (p)
            load_iota(p);
        Float eps = FloatLimits::min_invertible();
        for (int k = 0; k < n; k++) {
            if (p) { // Optionally pivot.
                int piv = find_pivot(a.diag()[Slice(k, n)]) + k;
                if (piv != k) {
                    std::swap(p[k], p[piv]);
                    a.swap_rows(k, piv);
                    a.swap_cols(k, piv);
                }
                if (k == 0)
                    eps = pre::abs(a(0, 0)) * FloatLimits::epsilon();
                // Positive semi-definite?
                if (!(pre::abs(a(k, k)) > eps)) {
                    for (int i = k; i < n; i++)
                        for (int j = i; j < n; j++)
                            a(i, j) = Field(0);
                    break;
                }
            }
            Field& akk = a(k, k);
            akk = pre::sqrt(akk);
            expect(pre::isfinite(akk) && pre::abs(akk) > eps,
                   "hermitian and positive semi-definite matrix");
            a(k, Slice(k + 1, n)) /= akk;
            for (int j = k + 1; j < n; j++)
                for (int i = k + 1; i < j + 1; i++) {
                    a(i, j) = a(i, j) - a(k, j) * pre::conj(a(k, i));
                    a(j, i) = pre::conj(a(i, j));
                }
        }
        force_upper_triangular(a);
    }

    /// Use factorization to solve linear system \f$ AX = B \f$.
    ///
    /// \param[in]  a  Matrix of coefficients, factorized in-place.
    /// \param[in]  b  Matrix of coefficients.
    /// \param[out] x  Matrix of solutions.
    /// \param[in]  p  _Optional_. Pivot vector.
    ///
    void chol_solve(
            MatView<const Field> a,
            MatView<const Field> b,
            MatView<Field> x,
            VecView<const int> p = {}) {
        ASSERT(a.is_square());
        ASSERT(b.rows() == a.rows() && //
               x.rows() == a.cols() && //
               x.cols() == b.cols());
        ASSERT(!p || a.size() == p.size());
        int m = b.rows();
        int n = b.cols();
        auto _ = cache_.scoped_push();
        auto y = cache_.template vector<Field>(m);
        for (int j = 0; j < n; j++) {
            // Solve R*y = b.
            for (int i = 0; i < m; i++) {
                auto k = Slice(0, i);
                y[i] = b(p ? p[i] : i, j) - dot(pre::conj(*a(k, i)), *y[k]);
                if (a(i, i) != Float(0))
                    y[i] /= pre::conj(a(i, i));
            }
            // Solve Rx = y.
            for (int i = m - 1; i >= 0; i--) {
                auto k = Slice(i + 1, m);
                y[i] -= dot(*a(i, k), *y[k]);
                if (a(i, i) != Float(0))
                    y[i] /= a(i, i);
            }
            for (int i = 0; i < m; i++)
                x(p ? p[i] : i, j) = y[i];
        }
    }

    /** \} */

  public:
    /// \name Factorization: QR
    /** \{ */

    /// Factorize.
    ///
    /// \param[inout] a  Matrix to factorize \f$ A \to R \f$.
    /// \param[out]   q  Matrix \f$ Q \f$.
    ///
    void qr(MatView<Field> a, MatView<Field> q) {
        ASSERT(q.size() == a.rows());
        ASSERT(q.is_square());
        load_identity(q);
        int m = a.rows();
        int n = a.cols();
        for (int k = 0; k < pre::min(m, n); k++)
            householderl(k, k, a, q);
        do_adjoint(q);
    }

  private:
    /// Apply Householder reflection from left.
    static constexpr void householderl_reflection(
            VecView<const Field> w, MatView<Field> a) noexcept {
        if (!a.empty()) {
            ASSERT(w.size() == a.rows());
            int m = a.rows();
            int n = a.cols();
            for (int j = 0; j < n; j++) {
                auto c = a(Slice(0, m), j);
                c = *c - (Float(2) * dot(pre::conj(*w), *c)) * (*w);
            }
        }
    }

    /// Householder reflection step targeting (0, 0) from left in submatrix.
    VecView<Field> householderl00(MatView<Field> a) {
        auto m = SliceInt(a.rows());
        auto n = SliceInt(a.cols());
        auto w = cache_.template vector<Field>(m);
        auto c = a(0 | m, 0);
        w = *c;
        c[0] = -norm2(w) * sign(w[0]);
        w[0] -= c[0];
        c(1 | m) = Field(0);
        normalize2(w);
        householderl_reflection(w, a(0 | m, 1 | n));
        return w;
    }

    /// Householder reflection step targeting (s, t) from left.
    void householderl(int s, int t, MatView<Field> x, MatView<Field> y) {
        auto m = SliceInt(x.rows());
        auto n = SliceInt(x.cols());
        if (s >= 0 && s < m && //
            t >= 0 && t < n) {
            auto _ = cache_.scoped_push();
            auto w = householderl00(x(s | m, t | n));
            if (!y.empty())
                householderl_reflection(w, y(s | m, 0 | m)); // Accumulate
        }
    }

    /// Householder reflection step targeting (s, t) from right.
    void householderr(int s, int t, MatView<Field> x, MatView<Field> y) {
        householderl(t, s, x.transpose(), y.transpose());
    }

    /** \} */

  public:
    /// \name Factorization: SVD
    /** \{ */

    /// Factorize.
    ///
    /// \param[inout] x  Matrix \f$ X \to \Sigma \f$.
    /// \param[out]   u  _Optional_. Matrix \f$ U \f$.
    /// \param[out]   v  _Optional_. Matrix \f$ V \f$.
    ///
    void svd(MatView<Field> x, MatView<Field> u = {}, MatView<Field> v = {}) {
        int m = x.rows();
        int n = x.cols();
        if (m < n) {
            svd(x.transpose(), v.transpose(), u.transpose());
            return;
        }
        // First reduce to upper bidiagonal.
        load_identity(u);
        load_identity(v);
        ASSERT(u.empty() || (u.is_square() && u.size() == x.rows()));
        ASSERT(v.empty() || (v.is_square() && v.size() == x.cols()));
        for (int k = 0; k < std::min(m, n); k++) {
            householderl(k, k + 0, x, u);
            householderr(k, k + 1, x, v);
        }
        // Golub-Kahan loop.
        int itr = 0;
        while (1) {
            // Find target indices.
            Float thresh = pre::max(pre::abs(*x.diag(0))) *
                           (FloatLimits::epsilon() * Float(8));
            int s = 0;
            while (s < n - 1 && pre::abs(x.diag(1)[s]) < thresh)
                s++;
            int t = s + 1;
            while (t < n - 1 && pre::abs(x.diag(1)[t]) > thresh)
                t++;
            if (t == n)
                break;
            // Form Gram submatrix terms.
            Field y0 = s + 1 < t ? x.diag(1)[t - 2] : Field();
            Field y1 = x.diag(1)[t - 1];
            Field z0 = x.diag(0)[t - 1];
            Field z1 = x.diag(0)[t];
            Float g00 = std::norm(y0) + std::norm(z0);
            Float g11 = std::norm(y1) + std::norm(z1);
            Float off = std::norm(z0) + std::norm(y1);
            // Form coefficients of characteristic polynomial.
            Float b = -g00 - g11;
            Float c = g00 * g11 - off;
            // Eigenvalues.
            b /= -2;
            Float d = std::fdim(b * b, c);
            Float lam0 = b + std::copysign(std::sqrt(d), b);
            Float lam1 = c / lam0;
            ASSERT(std::isfinite(lam0));
            ASSERT(std::isfinite(lam1));
            // Do Givens rotations.
            Field f = x.diag(0)[s] * x.diag(0)[s];
            Field g = x.diag(0)[s] * x.diag(1)[s];
            f -= std::abs(lam0 - g11) < std::abs(lam1 - g11) ? lam0 : lam1;
            for (int k = s; k < t; k++) {
                givensr(k, k + 1, f, g, x, v);
                f = x.diag(0)[k];
                g = x.diag(-1)[k];
                givensl(k, k + 1, f, g, x, u);
                if (k != t - 1) {
                    f = x.diag(1)[k];
                    g = x.diag(2)[k];
                }
            }
            expect(itr++ < 4096, "Golub-Kahan convergence");
        }
        do_adjoint(u);
        do_adjoint(v);
        // Force diagonal and positive.
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    x(i, j) = Field(0);
                }
                else {
                    if (x(i, i) != Float(0) && !v.empty())
                        v[i] *= pre::sign(x(i, i));
                    x(i, i) = std::abs(x(i, i));
                }
            }
    }

  private:
    /// Apply Givens rotation to pair of rows, given cosine and sine.
    static constexpr void givensl_rotation(
            int s,
            int t,
            Float cos_beta,
            Field sin_beta,
            MatView<Field> x) noexcept {
        for (int j = 0; j < int(x.cols()); j++) {
            Field x0 = x(s, j);
            Field x1 = x(t, j);
            x(s, j) = x0 * cos_beta + x1 * sin_beta;
            x(t, j) = x1 * cos_beta - x0 * pre::conj(sin_beta);
        }
    }

    /// Apply Givens rotation to pair of rows.
    static void givensl(
            int s,
            int t,
            Field f,
            Field g,
            MatView<Field> x,
            MatView<Field> y = {}) noexcept {
        ASSERT(y.empty() || (y.is_square()));
        if (s >= 0 && s < int(x.rows()) && //
            t >= 0 && t < int(x.rows())) {
            Float cos_beta = 1;
            Field sin_beta = 0;
            if (g != Float(0)) {
                if (f == Float(0)) {
                    cos_beta = 0;
                    sin_beta = pre::sign(g);
                }
                else {
                    Float absf = pre::abs(f);
                    Float absg = pre::abs(g);
                    cos_beta = absf;
                    sin_beta = (f / absf) * pre::conj(g);
                    Float denom = pre::hypot(absf, absg);
                    cos_beta /= denom;
                    sin_beta /= denom;
                }
            }
            givensl_rotation(s, t, cos_beta, sin_beta, x);
            givensl_rotation(s, t, cos_beta, sin_beta, y); // Accumulate
        }
    }

    /// Apply Givens rotation to pair of columns.
    static void givensr(
            int s,
            int t,
            Field f,
            Field g,
            MatView<Field> x,
            MatView<Field> y = {}) noexcept {
        return givensl(s, t, f, g, x.transpose(), y.transpose());
    }

    /** \} */
};

} // namespace pre
