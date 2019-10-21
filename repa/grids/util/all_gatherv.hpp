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

#include <boost/mpi.hpp>

namespace repa {
namespace util {

template <typename T, typename Size_type>
void all_gatherv_displ(boost::mpi::communicator &comm,
                 const std::vector<T> &sendbuf,
                 const std::vector<Size_type> &prefix_per_proc,
                 std::vector<T> &recvbuf)
{
    // MPI expects integer recvcounts and displacements for MPI_Allgatherv.
    // Copy idx_t vector to int vector.
    std::vector<int> recvcount(comm.size()), displ(comm.size());
    for (int i = 0; i < comm.size(); ++i) {
        recvcount[i] = static_cast<int>(prefix_per_proc[i + 1]
                                        - prefix_per_proc[i]);
        displ[i] = static_cast<int>(prefix_per_proc[i]);
    }

    MPI_Datatype mpi_type = boost::mpi::get_mpi_datatype(T{});

    MPI_Allgatherv(sendbuf.data(), static_cast<int>(sendbuf.size()), mpi_type,
                   recvbuf.data(), recvcount.data(), displ.data(), mpi_type,
                   comm);
}

} // namespace util
} // namespace repa
