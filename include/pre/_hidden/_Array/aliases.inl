/*-*- C++ -*-*/
#pragma once

namespace pre {

template <typename T>
using VecView = ArrayView<T, 1>;

template <typename T>
using Vec1 = Array<T, 1>;

template <typename T>
using Vec2 = Array<T, 2>;

template <typename T>
using Vec3 = Array<T, 3>;

template <typename T>
using Vec4 = Array<T, 4>;

template <typename T>
using MatView = ArrayView<T, 2>;

template <typename T>
using Mat1 = Array<T, 1, 1>;

template <typename T>
using Mat2 = Array<T, 2, 2>;

template <typename T>
using Mat3 = Array<T, 3, 3>;

template <typename T>
using Mat4 = Array<T, 4, 4>;

template <typename T>
using Mat1x1 = Array<T, 1, 1>;

template <typename T>
using Mat1x2 = Array<T, 1, 2>;

template <typename T>
using Mat1x3 = Array<T, 1, 3>;

template <typename T>
using Mat1x4 = Array<T, 1, 4>;

template <typename T>
using Mat2x1 = Array<T, 2, 1>;

template <typename T>
using Mat2x2 = Array<T, 2, 2>;

template <typename T>
using Mat2x3 = Array<T, 2, 3>;

template <typename T>
using Mat2x4 = Array<T, 2, 4>;

template <typename T>
using Mat3x1 = Array<T, 3, 1>;

template <typename T>
using Mat3x2 = Array<T, 3, 2>;

template <typename T>
using Mat3x3 = Array<T, 3, 3>;

template <typename T>
using Mat3x4 = Array<T, 3, 4>;

template <typename T>
using Mat4x1 = Array<T, 4, 1>;

template <typename T>
using Mat4x2 = Array<T, 4, 2>;

template <typename T>
using Mat4x3 = Array<T, 4, 3>;

template <typename T>
using Mat4x4 = Array<T, 4, 4>;

} // namespace pre
