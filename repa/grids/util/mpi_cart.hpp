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
#include "pargrid.hpp" // rank_type
#include <mpi.h>

namespace repa {
namespace util {

inline Vec3i mpi_cart_get_dims(MPI_Comm comm)
{
    Vec3i dims, _dummy, _dummy2;
    MPI_Cart_get(comm, 3, dims.data(), _dummy.data(), _dummy2.data());
    return dims;
}

inline Vec3i mpi_cart_get_coords(MPI_Comm comm)
{
    Vec3i coords, _dummy, _dummy2;
    MPI_Cart_get(comm, 3, _dummy.data(), _dummy2.data(), coords.data());
    return coords;
}

inline Vec3i mpi_cart_get_coords(MPI_Comm comm, rank_type r)
{
    Vec3i coords;
    MPI_Cart_coords(comm, static_cast<int>(r), 3, coords.data());
    return coords;
}

inline rank_type mpi_cart_rank(MPI_Comm comm, const Vec3i &coord)
{
    int r;
    MPI_Cart_rank(comm, coord.data(), &r);
    return static_cast<rank_type>(r);
}

inline int mpi_cart_get_color(MPI_Comm comm)
{
    Vec3i coords = mpi_cart_get_coords(comm);
    Vec3i dims = mpi_cart_get_dims(comm);

    for (int d = 0; d < 3; d++) {
        // Check if number of processes in this dimension is odd.
        // If this is the case, the nodes at the border of this dimension
        // must not be shifted.
        if ((dims[d] % 2 == 1) && (coords[d] == dims[d] - 1)) {
            return -1;
        }
    }

    return coords[0] % 2 + (coords[1] % 2) * 2 + (coords[2] % 2) * 4;
}

} // namespace util
} // namespace repa
