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

static MPI_Comm
make_cartesian_communicator(const boost::mpi::communicator &comm)
{
    Vec3i dims{0, 0, 0}, periods{1, 1, 1};
    MPI_Dims_create(comm.size(), 3, dims.data());

    MPI_Comm _comm_cart;
    MPI_Cart_create(comm, 3, dims.data(), periods.data(), true, &_comm_cart);

    return _comm_cart;
}

ParallelLCGrid::ParallelLCGrid(const boost::mpi::communicator &_comm,
                               Vec3d box_size,
                               double min_cell_size)
    : comm(_comm, boost::mpi::comm_duplicate),
      comm_cart(make_cartesian_communicator(_comm),
                boost::mpi::comm_take_ownership),
      box_l(box_size),
      node_grid(util::mpi_cart_get_dims(comm_cart)),
      node_pos(util::mpi_cart_get_coords(comm_cart)),
      max_range(min_cell_size)
{
}

} // namespace grids
} // namespace repa