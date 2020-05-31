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

#include "pargrid.hpp"

#include "grids/util/mpi_cart.hpp"

namespace repa {
namespace grids {

ParallelLCGrid::ParallelLCGrid(const boost::mpi::communicator &_comm,
                               Vec3d box_size,
                               double min_cell_size)
    : comm(_comm, boost::mpi::comm_duplicate),
      comm_cart(util::make_cartesian_communicator(_comm),
                boost::mpi::comm_take_ownership),
      box_size(box_size),
      min_cell_size(min_cell_size)
{
    // Grid should have at least three cells in each direction to be guaranteed
    // to work. Otherwise cell neighborhoods might be incorrect.
    //
    // Cell size is guaranteed to be:
    // min_cell_size <= cell_size <= 2 * min_cell_size
    ensure(box_size[0] >= 6 * min_cell_size && box_size[1] >= 6 * min_cell_size
               && box_size[2] >= 6 * min_cell_size,
           "Grid needs a minimum of three cells per dimension.");
}

} // namespace grids
} // namespace repa