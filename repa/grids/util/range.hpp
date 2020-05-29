/**
 * Copyright 2017-2020 Steffen Hirschmann
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

#include <boost/range/adaptors.hpp>
#include <boost/range/irange.hpp>

#include "strong_alias.hpp"

namespace repa {
namespace util {

/*
template <typename Int>
auto range(Int val) {
    return boost::irange(val);
}
*/

template <typename Int, typename Tag, Int Min, Int Max>
auto range(StrongAlias<Int, Tag, Min, Max> val)
{
    return boost::irange(static_cast<Int>(val))
           | boost::adaptors::transformed(
               [](Int i) { return StrongAlias<Int, Tag, Min, Max>{i}; });
}

} // namespace util
} // namespace repa