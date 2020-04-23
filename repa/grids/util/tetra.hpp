
/**
 * Copyright 2017-2020 Steffen Hirschmann
 * Copyright 2020 Benjamin Vier
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
#include <memory>

namespace repa {
namespace util {
namespace tetra {

// Opaque struct to reduce compile times on inclusion site
struct _Octagon_Impl;

extern int16_t precision;

// Initialize precision and box_size
void init_tetra(double min_cell_size, Vec3d box_size);

struct Octagon {
    Octagon();
    Octagon(const std::array<Vec3d, 8> &vertices, double max_cutoff = 0.0);
    Octagon(const Octagon &o) = delete;
    Octagon(Octagon &&o);
    ~Octagon();
    void operator=(Octagon o);

    /** Returns if this Octagon is valid and can be handled internally by this
     * module. Can only be called if a "max_cutoff" was passed to the
     * destructor. Otherwise, will throw a runtime_error.
     */
    bool is_valid() const;
    bool contains(const Vec3d &p) const;

private:
    std::unique_ptr<_Octagon_Impl> oi;

    friend void swap(Octagon &, Octagon &);
};

} // namespace tetra
} // namespace util
} // namespace repa