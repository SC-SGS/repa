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
#include <cmath>
#include <type_traits>

namespace repa {
namespace util {

template <typename T,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
T norm2(const T *v)
{
    return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

template <typename T,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
T norm2(const Vec3<T> &v)
{
    return norm2(v.data());
}

template <typename T,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
T dist2(const Vec3<T> &v, const Vec3<T> &w)
{
    Vec3<T> vw;
    for (typename Vec3<T>::size_type d = 0; d < v.size(); ++d)
        vw[d] = v[d] - w[d];
    return norm2(vw);
}

} // namespace util
} // namespace repa
