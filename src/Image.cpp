#include <iostream>
#include <pre-graphics/Image>
#include <pre-graphics/ImageSampler>

#define STB_IMAGE_IMPLEMENTATION 1
#define STB_IMAGE_WRITE_IMPLEMENTATION 1
#include "stb_image.h"
#include "stb_image_write.h"

namespace pre {

static void validate_axis(int axis) {
    if (axis < 0 or //
        axis > 3)
        throw std::invalid_argument("invalid axis");
}

std::vector<Image::LoadFunction> Image::load_functions;

std::vector<Image::SaveFunction> Image::save_functions;

void Image::resize(Dims new_dims, Type new_type, const void* new_data) {
    if (new_dims.prod() <= 0 or new_type == None) {
        clear();
        return;
    }
    type_ = new_type;
    dims_ = new_dims;
    data_.resize(dims_.prod() * Type_size(type_));
    if (new_data)
        std::copy(
                static_cast<const std::byte*>(new_data),
                static_cast<const std::byte*>(new_data) + data_.size(),
                data_.data());
}

void Image::resize(Dims new_dims) {
    if (new_dims.prod() <= 0) {
        clear();
        return;
    }
    Image new_image(new_dims, type_, nullptr, wrap);
    auto new_view = new_image.view();
    auto old_view = view();
    for (auto index : ArrayRange(ArrayIndex(min(dims_, new_dims)).join(1))) {
        auto* new_data = &new_view[index];
        auto* old_data = &old_view[index];
        std::copy(
                old_data, old_data + Type_size(type_), //
                new_data);
    }
    new_image.swap(*this);
}

template <concepts::image_data Data, size_t Channels>
static void resample1(
        ImageSampler1<Data, Channels> sampler0, //
        ImageSampler1<Data, Channels> sampler1) {
    int size0 = sampler0.view.size();
    int size1 = sampler1.view.size();
    if (size0 == size1) {
        sampler1.view = *sampler0.view;
        return;
    }
    float ratio = size0 / float(size1);
    float width = 2.0f * std::fmax(1.0f, ratio);
    for (int pos1 = 0; pos1 < size1; pos1++)
        sampler1.store(
                {pos1}, sampler0.filter(
                                {pos1 * ratio}, {width},
                                image_filter::mitchell<float>));
}

template <concepts::image_data Data, size_t Channels>
static void resample1(Image& image0, Image& image1, int axis, int& first) {
    auto sampler0 = make_sampler<Data, Channels>(image0, first);
    auto sampler1 = make_sampler<Data, Channels>(image1, first);
    sampler0.view = sampler0.view.transpose(2, axis);
    sampler1.view = sampler1.view.transpose(2, axis);
    std::swap(sampler0.wrap[2], sampler0.wrap[axis]);
    std::swap(sampler1.wrap[2], sampler1.wrap[axis]);
    first += Channels;
    int other0 = sampler0.view.sizes[0];
    int other1 = sampler0.view.sizes[1];
    for (int index0 = 0; index0 < other0; index0++)
        for (int index1 = 0; index1 < other1; index1++)
            resample1(
                    sampler0[index0][index1], //
                    sampler1[index0][index1]);
}

template <concepts::image_data Data>
static void resample1(Image& image0, int axis, int size) {
    Image::Type type = image0.type();
    Image::Dims dims = image0.dims();
    if (dims[axis] != size) {
        dims[axis] = size;
        Image image1(dims, type, nullptr, image0.wrap);
        for (int first = 0; first < dims[3];) {
            if (first + 4 <= dims[3])
                resample1<Data, 4>(image0, image1, axis, first);
            else if (first + 2 <= dims[3])
                resample1<Data, 2>(image0, image1, axis, first);
            else
                resample1<Data, 1>(image0, image1, axis, first);
        }
        image0.swap(image1);
    }
}

void Image::resample(Array<int, 3> new_dims) {
    if (empty() or new_dims.prod() <= 0) {
        clear();
        return;
    }
    Type prev_type = type_;
    packed_cast(Type_best_float(type_));
    for (int axis = 0; axis < 3; axis++) {
        switch (type_) {
        default:
        case Float32: resample1<float>(*this, axis, new_dims[axis]); break;
        case Float64: resample1<double>(*this, axis, new_dims[axis]); break;
        }
    }
    packed_cast(prev_type);
}

void Image::transpose(int axis0, int axis1) {
    validate_axis(axis0);
    validate_axis(axis1);
    if (axis0 == axis1)
        return;

    Dims new_dims = dims_;
    std::swap(new_dims[axis0], new_dims[axis1]);
    Image new_image(new_dims, type_, nullptr, wrap);
    new_image.view() = *view().transpose(axis0, axis1);
    new_image.swap(*this);
}

void Image::flip(int axis) {
    validate_axis(axis);
    ArrayRange<5> range;
    range.limit[0] = dims_[0];
    range.limit[1] = dims_[1];
    range.limit[2] = dims_[2];
    range.limit[3] = dims_[3];
    range.limit[4] = 1;
    range.limit[axis] = 1;
    for (auto index : range) {
        auto index0 = index;
        auto index1 = index;
        ssize_t& pos0 = index0[axis];
        ssize_t& pos1 = index1[axis];
        pos0 = 0;
        pos1 = dims_[axis] - 1;
        while (pos0 < pos1) {
            std::byte* ptr0 = &view()[index0];
            std::byte* ptr1 = &view()[index1];
            std::swap_ranges(ptr0, ptr0 + Type_size(type_), ptr1);
            ++pos0;
            --pos1;
        }
    }
}

void Image::swizzle(IteratorRange<const int*> ints) {
    if (empty() or ints.empty()) {
        clear();
        return;
    }
    for (int k = 0; k < int(ints.size()); k++)
        if (ints[k] < 0 or //
            ints[k] > dims_[3] - 1)
            throw std::invalid_argument("swizzle index out of range");
    Dims new_dims;
    new_dims[0] = dims_[0];
    new_dims[1] = dims_[1];
    new_dims[2] = dims_[2];
    new_dims[3] = ints.size();
    Image new_image(new_dims, type_, nullptr, wrap);
    auto go = [&]<typename Data>(Data) {
        auto new_view = new_image.view<Data>();
        auto old_view = view<Data>();
        for (auto index : ArrayRange<3>(new_dims)) {
            auto new_subview = new_view[index];
            auto old_subview = old_view[index];
            for (int k = 0; k < int(ints.size()); k++)
                new_subview[k] = old_subview[ints[k]];
        }
    };
    switch (Type_size(type_)) {
    default: break;
    case 1: go(std::uint8_t()); break;
    case 2: go(std::uint16_t()); break;
    case 4: go(std::uint32_t()); break;
    case 8: go(std::uint64_t()); break;
    }
    new_image.swap(*this);
}

void Image::stitch(int axis, const Image& other) {
    validate_axis(axis);
    auto& image0 = *this;
    auto& image1 = other;
    if (image0.empty()) {
        image0 = image1;
        return;
    }
    if (not(image0.type() == image1.type()))
        throw std::invalid_argument("incompatible image types");
    if (not(image0.dims() == image1.dims() or Array{0, 1, 2, 3} == axis).all())
        throw std::invalid_argument("incompatible image dimensions");

    Dims new_dims;
    new_dims = image0.dims();
    new_dims[axis] += image1.dims()[axis];
    Image new_image(new_dims, image0.type(), nullptr, image0.wrap);
    auto new_view = new_image.view().transpose(0, axis);
    auto image0_view = image0.view().transpose(0, axis);
    auto image1_view = image1.view().transpose(0, axis);
    int p0 = 0;
    int p1 = image0.dims()[axis];
    int p2 = image0.dims()[axis] + image1.dims()[axis];
    new_view(Slice(p0, p1)) = *image0_view;
    new_view(Slice(p1, p2)) = *image1_view;
    new_image.swap(*this);
}

void Image::srgb_encode(IteratorRange<const int*> ints) {
    if (empty() or ints.empty())
        return;
    Type prev_type = type_;
    packed_cast(Type_best_float(type_));
    auto go = [&]<typename Float>(Float) {
        for (auto i : ints)
            view<Float>().transpose(0, 3)[i].for_each([](Float& x) {
                x = pre::srgb_encode(pre::clamp(x, Float(0), Float(1)));
            });
    };
    switch (type_) {
    default:
    case Float32: go(float()); break;
    case Float64: go(double()); break;
    }
    packed_cast(prev_type);
}

void Image::srgb_decode(IteratorRange<const int*> ints) {
    if (empty() or ints.empty())
        return;
    Type prev_type = type_;
    packed_cast(Type_best_float(type_));
    auto go = [&]<typename Float>(Float) {
        for (auto i : ints)
            view<Float>().transpose(0, 3)[i].for_each([](Float& x) {
                x = pre::srgb_decode(pre::clamp(x, Float(0), Float(1)));
            });
    };
    switch (type_) {
    default:
    case Float32: go(float()); break;
    case Float64: go(double()); break;
    }
    packed_cast(prev_type);
}

template <Image::Type Source, Image::Type Target>
static constexpr void dispatch_cast(
        size_t num, const void* source, void* target) noexcept {
    using S = Image::astype<Source> const;
    using T = Image::astype<Target>;
    S* source_itr = static_cast<S*>(source);
    T* target_itr = static_cast<T*>(target);
    while (num-- != 0)
        if constexpr (
                (Source == Image::Float16 and (Target & Image::IsIntegral)) or
                (Target == Image::Float16 and (Source & Image::IsIntegral)))
            *target_itr++ = float(*source_itr++);
        else
            *target_itr++ = *source_itr++;
}

template <Image::Type Target>
static constexpr void dispatch_cast(
        size_t num,
        const Image::Type Source,
        const void* source,
        void* target) noexcept {
    switch (Source) {
#define SourceCase(T)                                                         \
    case T: dispatch_cast<T, Target>(num, source, target); break
        SourceCase(Image::Sint8);
        SourceCase(Image::Sint16);
        SourceCase(Image::Sint32);
        SourceCase(Image::Sint64);
        SourceCase(Image::Uint8);
        SourceCase(Image::Uint16);
        SourceCase(Image::Uint32);
        SourceCase(Image::Uint64);
        SourceCase(Image::Float16);
        SourceCase(Image::Float32);
        SourceCase(Image::Float64);
#undef SourceCase
    default: break;
    }
}

static constexpr void dispatch_cast(
        size_t num,
        const Image::Type Source,
        const Image::Type Target,
        const void* source,
        void* target) noexcept {
    switch (Target) {
#define TargetCase(T)                                                         \
    case T: dispatch_cast<T>(num, Source, source, target); break
        TargetCase(Image::Sint8);
        TargetCase(Image::Sint16);
        TargetCase(Image::Sint32);
        TargetCase(Image::Sint64);
        TargetCase(Image::Uint8);
        TargetCase(Image::Uint16);
        TargetCase(Image::Uint32);
        TargetCase(Image::Uint64);
        TargetCase(Image::Float16);
        TargetCase(Image::Float32);
        TargetCase(Image::Float64);
#undef TargetCase
    default: break;
    }
}

void Image::cast(Type new_type) {
    Type cur_type = type_;
    if (cur_type == new_type or //
        cur_type == None or new_type == None)
        return;
    if ((cur_type & IsIntegral) and //
        (new_type & IsIntegral) and //
        Type_size(cur_type) == Type_size(new_type)) {
        type_ = new_type;
        return;
    }
    Image tmp = std::move(*this);
    resize(tmp.dims(), new_type);
    dispatch_cast(
            dims().prod(), //
            cur_type, new_type, tmp.data(), data());
}

void Image::packed_cast(Type new_type) {
    Type cur_type = type_;
    if (cur_type == new_type or //
        cur_type == None or new_type == None)
        return;
    // Casting between integral types?
    if ((cur_type & IsIntegral) and //
        (new_type & IsIntegral)) {
        cast(Float64);
        packed_cast(new_type);
        return;
    }
    // Casting between floating types?
    if (not(cur_type & IsIntegral) and //
        not(new_type & IsIntegral)) {
        cast(new_type); // Just cast
        return;
    }
    // Casting floating to integral?
    if (not(cur_type & IsIntegral)) {
        if (cur_type == Float16) {
            cur_type = Float32;
            cast(Float32); // If half, cast up to float.
        }
        double vfac = 0;
        switch (new_type) {
        default:
        case Sint8: vfac = Maximum<std::int8_t>; break;
        case Sint16: vfac = Maximum<std::int16_t>; break;
        case Sint32: vfac = Maximum<std::int32_t>; break;
        case Sint64: vfac = Maximum<std::int64_t>; break;
        case Uint8: vfac = Maximum<std::uint8_t>; break;
        case Uint16: vfac = Maximum<std::uint16_t>; break;
        case Uint32: vfac = Maximum<std::uint32_t>; break;
        case Uint64: vfac = Maximum<std::uint64_t>; break;
        }
        if (cur_type == Float32) {
            float vmin = (new_type & IsUnsigned) ? 0 : -1;
            float vmax = 1;
            for (float& v : IteratorRange(
                         data<float>(), //
                         data<float>() + dims().prod())) {
                v = std::max(v, vmin);
                v = std::min(v, vmax);
                v *= vfac;
            }
        }
        else {
            double vmin = (new_type & IsUnsigned) ? 0 : -1;
            double vmax = 1;
            for (double& v : IteratorRange(
                         data<double>(), //
                         data<double>() + dims().prod())) {
                v = std::max(v, vmin);
                v = std::min(v, vmax);
                v *= vfac;
            }
        }
        cast(new_type);
    }
    // Casting integral to floating.
    else {
        double vfac = 0;
        switch (cur_type) {
        default:
        case Sint8: vfac = 1.0 / Maximum<std::int8_t>; break;
        case Sint16: vfac = 1.0 / Maximum<std::int16_t>; break;
        case Sint32: vfac = 1.0 / Maximum<std::int32_t>; break;
        case Sint64: vfac = 1.0 / Maximum<std::int64_t>; break;
        case Uint8: vfac = 1.0 / Maximum<std::uint8_t>; break;
        case Uint16: vfac = 1.0 / Maximum<std::uint16_t>; break;
        case Uint32: vfac = 1.0 / Maximum<std::uint32_t>; break;
        case Uint64: vfac = 1.0 / Maximum<std::uint64_t>; break;
        }
        if (Type_size(cur_type) < sizeof(float)) {
            cast(Float32);
            for (float& v : IteratorRange(
                         data<float>(), //
                         data<float>() + dims().prod()))
                v *= vfac;
        }
        else {
            cast(Float64);
            for (double& v : IteratorRange(
                         data<double>(), //
                         data<double>() + dims().prod()))
                v *= vfac;
        }
        cast(new_type);
    }
}

bool Image::load(const std::string& filename) try {
    clear();
    for (const LoadFunction& func : load_functions)
        if (func(filename, *this))
            return true;
    return false;
}
catch (...) {
    std::cerr << "Error in Image::load()!\n";
    throw;
}

bool Image::save(const std::string& filename) const try {
    for (const SaveFunction& func : save_functions)
        if (func(filename, *this))
            return true;
    return false;
}
catch (...) {
    std::cerr << "Error in Image::save()!\n";
    throw;
}

static bool Image_stbi_load(const std::string& filename, Image& im) {
    void* data = nullptr;
    int sizex = 0;
    int sizey = 0;
    int bands = 0;
    Image::Type type = Image::None;
    if (stbi_is_hdr(filename.c_str())) {
        data = stbi_loadf(filename.c_str(), &sizex, &sizey, &bands, 0);
        type = Image::Float32;
    }
    else if (stbi_is_16_bit(filename.c_str())) {
        data = stbi_load_16(filename.c_str(), &sizex, &sizey, &bands, 0);
        type = Image::Uint16;
    }
    else {
        data = stbi_load(filename.c_str(), &sizex, &sizey, &bands, 0);
        type = Image::Uint8;
    }
    if (data) { // Success?
        im.resize({1, sizey, sizex, bands}, type, data);
        stbi_image_free(data);
        return true;
    }
    else {
        return false;
    }
}

static bool Image_stbi_save(const std::string& filename, const Image& im) {
    if (im.type() != Image::Uint8)
        return false;
    Image::Dims dims = im.dims();
    if (dims[0] != 1)
        return false;
    int sizey = dims[1];
    int sizex = dims[2];
    int bands = dims[3];
    if (filename.ends_with(".png"))
        return stbi_write_png(
                filename.c_str(), sizex, sizey, bands, im.data(),
                sizex * bands);
    if (filename.ends_with(".bmp"))
        return stbi_write_bmp(
                filename.c_str(), sizex, sizey, bands, im.data());
    if (filename.ends_with(".tga"))
        return stbi_write_tga(
                filename.c_str(), sizex, sizey, bands, im.data());
    if (filename.ends_with(".jpg") or filename.ends_with(".jpeg"))
        return stbi_write_jpg(
                filename.c_str(), sizex, sizey, bands, im.data(), 90);
    return false;
}

struct Stbi {
    Stbi() {
        Image::register_load(Image_stbi_load);
        Image::register_save(Image_stbi_save);
    }
};

static Stbi stbi_init;

} // namespace pre
