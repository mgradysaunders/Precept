#include <pre-graphics/Converger>
#include <pre-graphics/Microsurface>

namespace pre {

namespace microsurface {

Vec2<double> TrowbridgeReitzSlope::P11_sample(
        double u0, double u1, double cos_thetao) const noexcept {
    if (cos_thetao > 0.99999) {
        double r = pre::sqrt(u0 / (1 - u0));
        double phi = 2 * pi * u1;
        return {r * pre::cos(phi), //
                r * pre::sin(phi)};
    }
    else {
        double sin_thetao = pre::sqrt(1 - cos_thetao * cos_thetao);
        double tan_thetao = sin_thetao / cos_thetao;
        double mu = u0 * (1 + 1 / cos_thetao) - 1;
        double nu = 1 / (1 - mu * mu);
        double q = pre::sqrt(pre::fmax(
                0.0, mu * mu * nu - nu * (1 - nu) * tan_thetao * tan_thetao));
        double t0 = -nu * tan_thetao - q;
        double t1 = -nu * tan_thetao + q;
        double m0 = mu < 0 or t1 * sin_thetao > cos_thetao ? t0 : t1;
        double m1 = 1;
        if (u1 > 0.5) {
            u1 = 2 * u1 - 1;
        }
        else {
            u1 = 1 - 2 * u1;
            m1 = -1;
        }
        m1 *= pre::sqrt(1 + m0 * m0) *
              (u1 * (u1 * (u1 * 0.273850 - 0.733690) + 0.463410)) /
              (u1 * (u1 * (u1 * 0.093073 + 0.309420) - 1.000000) + 0.597999);
        return {m0, m1};
    }
}

Vec2<double> BeckmannSlope::P11_sample(
        double u0, double u1, double cos_thetao) const noexcept {
    if (cos_thetao > 0.99999) {
        double r = std::sqrt(-std::log1p(-u0));
        double phi = 2 * pi * u1;
        return {r * std::cos(phi), //
                r * std::sin(phi)};
    }
    else {
        u0 = std::fmax(u0, 1e-6);
        double sin_thetao = pre::sqrt(1 - cos_thetao * cos_thetao);
        double cot_thetao = cos_thetao / sin_thetao;
        auto c11 = [=](double a) {
            return 0.5 * inv_sqrtpi * sin_thetao * std::exp(-a * a) +
                   0.5 * cos_thetao * std::erfc(-a);
        };
        double cnorm = 1 / c11(cot_thetao);
        if (cnorm > 1e+6 or not std::isfinite(cnorm))
            return {};
        double xmin = -1;
        double xmax = std::erf(cot_thetao);
        Converger<double> converger;
        converger.max_iters = 10;
        converger.lower_bound = xmin;
        converger.upper_bound = xmax;
        converger.target = u0;
        converger.cutoff = 1e-6;
        auto f = [&](double x) {
            double a = erfinv(x);
            return x >= cot_thetao ? 1 : cnorm * c11(a);
        };
        auto g = [&](double x) {
            double a = erfinv(x);
            return cnorm * (cos_thetao - a * sin_thetao) / 2;
        };
        // Improved initial guess lifted from PBRT-v3 source.
        double thetao = std::acos(cos_thetao);
        double x = -0.0564;
        x = std::fma(thetao, x, +0.4265);
        x = std::fma(thetao, x, -0.876);
        x = std::fma(thetao, x, +1.0);
        x = xmax - (1 + xmax) * std::pow(1 - u0, x);
        // Do numerical inversion.
        if (converger(x, f, g))
            return {erfinv(x), erfinv(2 * u1 - 1)};
    }
    return {};
}

UniformHeight uniform_height;

NormalHeight normal_height;

TrowbridgeReitzSlope trowbridge_reitz_slope;

BeckmannSlope beckmann_slope;

Fresnel fresnel;

} // namespace microsurface

Microsurface::Result Microsurface::simulate_path(
        const Vec3<double>& wo,
        const Vec3<double>& wi,
        int min_order,
        int max_order,
        bool importance_mode) const noexcept {
    Result result;
    double energyk = 1;
    double heightk = 1 + height->C1inv(0.99999);
    Vec3<double> wk = -wo;
    bool wo_outside = wo[2] > 0;
    bool wi_outside = wi[2] > 0;
    bool wk_outside = wo_outside;
    if (!wk_outside)
        heightk = -heightk; // Flip height
    double p0 = 0;
    for (int order = 0; max_order <= 0 or order < max_order; order++) {
        float u = generate_canonical<float>(random);
        heightk = wk_outside //
                          ? +h_sample(u, +wk, +heightk)
                          : -h_sample(u, -wk, -heightk);
        if (std::isinf(heightk))
            break;
        if (order >= min_order) {
            double ek = energyk;
            double pk =
                    phase(-wk, wi,    //
                          wk_outside, //
                          wi_outside, //
                          &ek, importance_mode);
            double fk = pk;
            fk *= wi_outside ? G1(+wi, +heightk) : G1(-wi, -heightk);
            fk *= order == 0 ? 0.5 : p0 / (p0 + pk); // MIS weight.
            ek *= fk;
            if (std::isfinite(fk)) {
                result.f += ek;
                result.fpdf += fk;
            }
        }
        wk = phase_sample(
                -wk, wk_outside, &wk_outside, &energyk, importance_mode);
        // Remember MIS weight term.
        if (order == 0)
            p0 = phase(
                    wk, wo, wk_outside, wo_outside, nullptr, importance_mode);
        // Safety check.
        if (not pre::isfinite(heightk) or //
            not pre::isfinite(energyk) or energyk == 0 or wk[2] == 0)
            break;
    }
    return result;
}

double LambertianMicrosurface::phase(
        const Vec3<double>& wo,
        const Vec3<double>& wi,
        bool wo_outside,
        bool wi_outside,
        double* energy,
        bool) const noexcept {
    using microsurface::inv_pi;
    float u0 = generate_canonical<float>(random);
    float u1 = generate_canonical<float>(random);
    Vec3<double> wm = wo_outside //
                              ? +Dwo_sample(u0, u1, +wo)
                              : -Dwo_sample(u0, u1, -wo);
    if (energy)
        *energy *= r0 + t0;
    return wo_outside == wi_outside //
                   ? inv_pi * pre::max(dot(+wm, wi), 0.0) * r0 / (r0 + t0)
                   : inv_pi * pre::max(dot(-wm, wi), 0.0) * t0 / (r0 + t0);
}

Vec3<double> LambertianMicrosurface::phase_sample(
        Vec3<double> wo,
        bool wo_outside,
        bool* wi_outside,
        double* energy,
        bool) const noexcept {
    float u0 = generate_canonical<float>(random);
    float u1 = generate_canonical<float>(random);
    Vec3<double> wm = wo_outside //
                              ? +Dwo_sample(u0, u1, +wo)
                              : -Dwo_sample(u0, u1, -wo);
    Vec3<double> wi = //
            Vec3<double>::cosine_hemisphere_pdf_sample(
                    {generate_canonical<float>(random),
                     generate_canonical<float>(random)});
    if (generate_canonical<float>(random) < r0 / (r0 + t0))
        *wi_outside = wo_outside;
    else
        *wi_outside = not wo_outside, wi[2] *= -1;
    *energy *= r0 + t0;
    return normalize(dot(Mat3<double>::build_onb(wm), wi));
}

double DielectricMicrosurface::phase(
        const Vec3<double>& wo,
        const Vec3<double>& wi,
        bool wo_outside,
        bool wi_outside,
        double* energy,
        bool importance_mode) const noexcept {
    double etao = wo_outside ? eta_above : eta_below;
    double etat = wo_outside ? eta_below : eta_above;
    double eta = etao / etat;
    if (wo_outside == wi_outside) {
        Vec3<double> wm = normalize(wo + wi);
        double cos_thetao = dot(wo, wm);
        auto [Fr, Ft] = fresnel->F(etao, etat, cos_thetao);
        Fr *= Fr0;
        Ft *= Ft0;
        if (energy)
            *energy *= Fr + Ft;
        double D = Dwo(                //
                wo_outside ? wo : -wo, //
                wo_outside ? wm : -wm);
        return finite_or_zero(D / (4 * cos_thetao) * Fr / (Fr + Ft));
    }
    else {
        Vec3<double> vm = eta * wo + wi;
        if (vm[2] < 0)
            vm = -vm;
        if (not wo_outside)
            vm = -vm;
        Vec3<double> wm = normalize(vm);
        double cos_thetao = dot(wo, wm);
        double cos_thetai = dot(wi, wm);
        if (not(cos_thetao > 0 and //
                cos_thetai < 0))
            return 0;
        auto [Fr, Ft] = fresnel->F(etao, etat, cos_thetao);
        Fr *= Fr0;
        Ft *= Ft0;
        if (energy) {
            *energy *= Fr + Ft;
            if (importance_mode)
                *energy *= eta * eta;
        }
        double D = Dwo(                //
                wo_outside ? wo : -wo, //
                wo_outside ? wm : -wm);
        return finite_or_zero(D * -cos_thetai / dot(vm, vm) * Ft / (Fr + Ft));
    }
}

Vec3<double> DielectricMicrosurface::phase_sample(
        Vec3<double> wo,
        bool wo_outside,
        bool* wi_outside,
        double* energy,
        bool importance_mode) const noexcept {
    double etao = wo_outside ? eta_above : eta_below;
    double etat = wo_outside ? eta_below : eta_above;
    double eta = etao / etat;
    float u0 = generate_canonical<float>(random);
    float u1 = generate_canonical<float>(random);
    Vec3<double> wm = wo_outside //
                              ? +Dwo_sample(u0, u1, +wo)
                              : -Dwo_sample(u0, u1, -wo);
    double cos_thetao = dot(wo, wm);
    double cos_thetat = 0;
    auto [Fr, Ft] = fresnel->F(etao, etat, cos_thetao, &cos_thetat);
    Fr *= Fr0;
    Ft *= Ft0;
    *energy *= Fr + Ft;
    if (generate_canonical<float>(random) < Fr / (Fr + Ft)) {
        *wi_outside = wo_outside;
        return normalize(2 * cos_thetao * wm - wo);
    }
    else {
        if (importance_mode)
            *energy *= eta * eta;
        *wi_outside = not wo_outside;
        return normalize((eta * cos_thetao + cos_thetat) * wm - eta * wo);
    }
}

} // namespace pre
