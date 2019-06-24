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

#include "common_types.hpp"
#include <array>

namespace repa {
namespace util {

struct NeighborOffsets3D {
    // Offsets of cell and process neighbors -- respects cell ordering
    // requirements of pargrid.hpp: 0 cell itself 1-13 hs neigh 14-26 fs neigh
    static const std::array<Vec3i, 27> raw;
};

} // namespace util
} // namespace repa
