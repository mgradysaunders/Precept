#include <pre-graphics/Converger>
#include <pre-graphics/AnimationCurve>

namespace pre {

float AnimationCurve::operator()(float time) const {
    switch (keyframes.size()) {
    case 0: return 0;
    case 1: return keyframes.front().value;
    default: break;
    }
    float time0 = keyframes.front().time;
    float time1 = keyframes.back().time;
    if (time < time0) {
        switch (wrap_mode_before) {
        default:
        case Wrap::Clamp: return keyframes.front().value;
        case Wrap::Repeat: time = frepeat(time, time0, time1); break;
        case Wrap::Mirror: time = fmirror(time, time0, time1); break;
        }
    }
    else if (time > time1) {
        switch (wrap_mode_after) {
        default:
        case Wrap::Clamp: return keyframes.back().value;
        case Wrap::Repeat: time = frepeat(time, time0, time1); break;
        case Wrap::Mirror: time = fmirror(time, time0, time1); break;
        }
    }
    const Keyframe* key0 = nullptr;
    const Keyframe* key1 = nullptr;
    if (keyframes.size() == 2) {
        key0 = &keyframes.front();
        key1 = &keyframes.back();
    }
    else {
        auto itr = std::lower_bound(
                keyframes.begin(), //
                keyframes.end(), time,
                [](const Keyframe& keyframe, float time) {
                    return keyframe.time < time;
                });
        if (itr == keyframes.begin())
            return keyframes.front().value;
        if (itr == keyframes.end())
            return keyframes.back().value;
        key0 = &*(itr - 1);
        key1 = &*(itr - 0);
    }
    time0 = key0->time;
    time1 = key1->time;
    float t = (time - time0) / (time1 - time0);
    float x0 = time0;
    float x1 = time1;
    float y0 = key0->value, m0 = key0->out_slope;
    float y1 = key1->value, m1 = key1->in_slope;
    if ((int(key0->weight_mode) & 2) == 0 and
        (int(key1->weight_mode) & 1) == 0)
        return bezier(
                t, y0,                   //
                y0 + m0 * (x1 - x0) / 3, //
                y1 - m1 * (x1 - x0) / 3, y1);

    // Weighted evaluation.
    float dx0 = !(int(key0->weight_mode) & 2) ? 0.333333f : key0->out_weight;
    float dx1 = !(int(key1->weight_mode) & 1) ? 0.333333f : key1->in_weight;
    dx0 *= x1 - x0;
    dx1 *= x1 - x0;
    float x = lerp(t, x0, x1);
    float xs[4] = {x0, x0 + dx0, x1 - dx1, x1};
    float ys[4] = {y0, y0 + m0 * dx0, y1 - m1 * dx1, y1};
    float ds[3] = {dx0, xs[2] - xs[1], dx1};

    // Solve with Newton-Raphson iteration.
    Converger<float> converger;
    converger.max_iters = 10;
    converger.lower_bound = 0;
    converger.upper_bound = 1;
    converger.target = x;
    converger.cutoff = 1e-4f;
    auto f = [&](float u) { return bezier(u, xs[0], xs[1], xs[2], xs[3]); };
    auto g = [&](float u) { return bezier(u, ds[0], ds[1], ds[2]) * 3; };
    if (not converger(t, f, g))
        throw std::runtime_error(
                "AnimationCurve::operator() failed to converge!");

    return bezier(t, ys[0], ys[1], ys[2], ys[3]);
}

} // namespace pre
