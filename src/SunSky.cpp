#include <pre-graphics/SunSky>

namespace pre {

Vec3<double> SunSky::sun_direction(
    double lat, double lon, double day) noexcept {
    double sin_latitude = std::sin(lat);
    double cos_latitude = std::cos(lat);
    double n = 360.0_degrees / 365.24;
    double d = std::floor(day) + 1;
    double a = n * (d + 9);
    double b = a + std::sin(n * (d - 3)) * 1.914_degrees;
    double c = a - std::atan(std::tan(b) / std::cos(23.44_degrees));
    c /= 180.0_degrees;
    c -= std::round(c);
    double t = c / 2; // Solar time
    t += fract(day);
    t -= fract(lon / 15.0_degrees) / 24 - 0.5 / 24;
    double omega = 15.0_degrees * (24 * t - 12); // Hour angle
    double sin_omega = std::sin(omega);
    double cos_omega = std::cos(omega);
    double sin_delta = std::cos(b) * -std::sin(23.44_degrees); // Declination
    double cos_delta = std::sqrt(std::fmax(0.0, 1.0 - sin_delta * sin_delta));
    return normalize(
        -cos_delta * sin_omega,
        cos_latitude * sin_delta - sin_latitude * cos_delta * cos_omega,
        sin_latitude * sin_delta + cos_latitude * cos_delta * cos_omega);
}

SunSky::DawnDusk SunSky::sun_dawn_dusk(
    double lat, double lon, int day) noexcept {
    double tan_latitude = std::tan(lat);
    double n = 360.0_degrees / 365.24;
    double d = day + 1;
    double a = n * (d + 9);
    double b = a + std::sin(n * (d - 3)) * 1.914_degrees;
    double c = a - std::atan(std::tan(b) / std::cos(23.44_degrees));
    c /= 180.0_degrees;
    c -= std::round(c);
    double sin_delta = std::cos(b) * -std::sin(23.44_degrees); // Declination
    double cos_delta = std::sqrt(std::fmax(0.0, 1.0 - sin_delta * sin_delta));
    double cos_omega = -tan_latitude * sin_delta / cos_delta;
    if (!(cos_omega > -0.9))
        return {0, 1};
    if (!(cos_omega < +0.9))
        return {0, 0};
    double omega = std::acos(cos_omega);
    double t0 = -omega / (15.0_degrees * 24) + 0.5;
    double t1 = +omega / (15.0_degrees * 24) + 0.5;
    t0 = fract(t0 + fract(lon / 15.0_degrees) / 24 - 0.5 / 24 - c / 2);
    t1 = fract(t1 + fract(lon / 15.0_degrees) / 24 - 0.5 / 24 - c / 2);
    if (t0 > t1)
        std::swap(t0, t1);
    return {t0, t1};
}

} // namespace pre
