
#pragma once

#include "common_types.hpp"
#include <type_traits>

namespace repa {
namespace util {

template <typename T>
Vec3<T> vadd(const Vec3<T> &a, const Vec3<T> &b)
{
    return {{a[0] + b[0], a[1] + b[1], a[2] + b[2]}};
}

template <typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
Vec3<T> vadd_mod(const Vec3<T> &a, const Vec3<T> &b, const Vec3<T> &mod)
{
    auto r = vadd(a, b);

    for (typename Vec3<T>::size_type i = 0; i < a.size(); ++i) {
        if (r[i] < T(0) || r[i] >= mod[i])
            r[i] -= (r[i] / mod[i]) * mod[i];
    }
    return r;
}

} // namespace util
} // namespace repa