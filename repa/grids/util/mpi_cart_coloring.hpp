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

#include "mpi_cart.hpp"

namespace repa {
namespace util {

namespace __impl {
/** Returns a color for independent set traveral.
 * All neighboring processes in a Cartesian setting will have different colors.
 * There are at most 8 colors.
 * To guarantee 8 colors, If the number of processes is odd in a direction,
 * returns -1 on all processes with the highest coordinate in this direction.
 */
inline int get_color(MPI_Comm comm)
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
} // namespace __impl

struct independent_process_sets {
    independent_process_sets(const boost::mpi::communicator &comm_cart)
        : col(__impl::get_color(comm_cart))
    {
        assert(comm_cart.has_cartesian_topology());
    }

    template <typename Func>
    independent_process_sets &for_each(Func &&f)
    {
        _round_func = f;
        return *this;
    }

    template <typename Func>
    independent_process_sets &for_all_after_each_round(Func &&f)
    {
        _after_func = f;
        return *this;
    }

    void operator()()
    {
        for (int color = 0; color < 8; color++) {
            if (color == col) {
                _round_func();
            }
            _after_func();
        }
    }

private:
    const int col;
    std::function<void(void)> _round_func, _after_func;
};

} // namespace util
} // namespace repa