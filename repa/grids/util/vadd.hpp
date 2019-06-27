/**
 * Copyright 2017-2019 Steffen Hirschmann
 *
 * This file is part of Repa.

 * Repa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * Repa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with Repa.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include "common_types.hpp"
#include <type_traits>
#include <cmath>

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
            r[i] -= static_cast<T>(std::floor(static_cast<double>(r[i]) / mod[i])) * mod[i];
    }
    return r;
}

} // namespace util
} // namespace repa