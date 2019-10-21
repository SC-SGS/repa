/**
 * Copyright 2017-2019 Steffen Hirschmann
 *
 * This file is part of Repa.
 *
 * Repa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Repa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Repa.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include <array>
#include <functional>
#include <type_traits>
#include <vector>

namespace repa {

template <typename T>
using Vec3 = std::array<T, 3>;

typedef Vec3<int> Vec3i;
typedef Vec3<double> Vec3d;

typedef std::function<std::vector<double>(void)> CellMetric;
typedef std::function<double(int, int)> CellCellMetric;
typedef std::function<void(void)> Thunk;

/** Type that behaves like an integral POD and can be restricted to a range.
 * Range is tested at construction time.
 * Should not compile to any overhead on reasonable compilers and optimization
 * levels.
 */
template <typename T,
          T min,
          T max,
          typename = typename std::enable_if_t<std::is_integral<T>::value>>
struct IntegralRange {
    typedef T value_type;

    template <typename S,
              typename
              = typename std::enable_if_t<std::is_integral<S>::value>>
    inline IntegralRange(S v) : value(static_cast<T>(v))
    {
        // Evaluate range check on wide base type.
        typedef std::common_type_t<T, S> base_type;
        assert(static_cast<base_type>(v) >= static_cast<base_type>(min)
               && static_cast<base_type>(v) <= static_cast<base_type>(max));
    }
    inline operator value_type()
    {
        return value;
    }

private:
    value_type value;
};

typedef IntegralRange<std::int_fast32_t, 0, 26> fs_neighidx;

} // namespace repa