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

#include "mpi_graph.hpp"
#include <boost/mpi.hpp>

namespace repa {
namespace util {

namespace __impl {

/** MPI_Scatter "data" to the processes indicated by "neighbors" as well as
 * MPI_Gather from them and return this data.
 * Corresponds to a MPI_Alltoall among "neighbors".
 */
template <typename T>
std::vector<T> mpi_subset_alltoall(const boost::mpi::communicator &comm,
                                   const std::vector<rank_type> &neighbors,
                                   const std::vector<T> &data)
{
    std::vector<boost::mpi::request> sreq_cells(neighbors.size());
    std::vector<boost::mpi::request> rreq_cells(neighbors.size());

    for (rank_index_type i = 0; i < neighbors.size(); ++i) {
        sreq_cells[i] = comm.isend(neighbors[i], 2, data[i]);
    }

    std::vector<T> gathered_data(neighbors.size());
    for (rank_index_type i = 0; i < neighbors.size(); ++i) {
        rreq_cells[i] = comm.irecv(neighbors[i], 2, gathered_data[i]);
    }

    boost::mpi::wait_all(std::begin(rreq_cells), std::end(rreq_cells));
    boost::mpi::wait_all(std::begin(sreq_cells), std::end(sreq_cells));

    return gathered_data;
}

} // namespace __impl

/** MPI_Neighbor_alltoall.
 */
template <typename T>
std::vector<T> mpi_neighbor_alltoall(const boost::mpi::communicator &neighcomm,
                                     const std::vector<T> &data)
{
    assert(has_dist_graph_topology(neighcomm));
    return __impl::mpi_subset_alltoall(
        neighcomm, mpi_undirected_neighbors(neighcomm), data);
}

} // namespace util
} // namespace repa
