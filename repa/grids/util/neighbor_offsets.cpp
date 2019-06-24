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

#include "neighbor_offsets.hpp"

namespace repa {
namespace util {
// clang-format off
const std::array<Vec3i, 27> NeighborOffsets3D::raw = {{
    // Cell itself
    {{0, 0, 0}},
    // Half shell begin
    {{1, 0, 0}},
    {{-1, 1, 0}},
    {{0, 1, 0}},
    {{1, 1, 0}},
    {{-1, -1, 1}},
    {{0, -1, 1}},
    {{1, -1, 1}},
    {{-1, 0, 1}},
    {{0, 0, 1}},
    {{1, 0, 1}},
    {{-1, 1, 1}},
    {{0, 1, 1}},
    {{1, 1, 1}},
    // Full shell begin
    {{-1, -1, -1}},
    {{0, -1, -1}},
    {{1, -1, -1}},
    {{-1, 0, -1}},
    {{0, 0, -1}},
    {{1, 0, -1}},
    {{-1, 1, -1}},
    {{0, 1, -1}},
    {{1, 1, -1}},
    {{-1, -1, 0}},
    {{0, -1, 0}},
    {{1, -1, 0}},
    {{-1, 0, 0}}}};
// clang-format on
} // namespace util
} // namespace repa
