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
#include <mpi.h>

namespace repa {
namespace util {

/** Returns the number of neighbors in an undirected Dist_graph communicator.
 * This function simply assumes that the graph topology is undirected and
 * does not verify it!
 */
inline int mpi_undirected_neighbor_count(MPI_Comm neighcomm)
{
    int indegree = 0, outdegree = 0, weighted = 0;
    MPI_Dist_graph_neighbors_count(neighcomm, &indegree, &outdegree, &weighted);
    return indegree;
}

/** Returns a boost::mpi::communicator with a undirected Dist_graph
 * topology.
 */
inline boost::mpi::communicator
directed_mpi_communicator(MPI_Comm parent_communicator,
                          const std::vector<int> &source_ranks,
                          const std::vector<int> &destination_ranks)
{
    MPI_Comm comm;
    // Edges to all processes in "neighbors"
    MPI_Dist_graph_create_adjacent(
        parent_communicator, source_ranks.size(), source_ranks.data(),
        static_cast<const int *>(MPI_UNWEIGHTED), destination_ranks.size(),
        destination_ranks.data(), static_cast<const int *>(MPI_UNWEIGHTED),
        MPI_INFO_NULL, 0, &comm);

#ifndef NDEBUG
    int indegree = 0, outdegree = 0, weighted = 0;
    MPI_Dist_graph_neighbors_count(comm, &indegree, &outdegree, &weighted);
    assert(!weighted);
    assert(static_cast<size_t>(indegree) == source_ranks.size());
    assert(static_cast<size_t>(outdegree) == destination_ranks.size());
    std::vector<int> __ineighs(indegree, -1), __iw(indegree, -1);
    std::vector<int> __oneighs(outdegree, -1), __ow(outdegree, -1);
    MPI_Dist_graph_neighbors(comm, indegree, __ineighs.data(), __iw.data(),
                             outdegree, __oneighs.data(), __ow.data());
    for (size_t i = 0; i < source_ranks.size(); ++i)
        assert(__ineighs[i] == source_ranks[i]);
    for (size_t i = 0; i < destination_ranks.size(); ++i)
        assert(__oneighs[i] == destination_ranks[i]);
#endif

    return boost::mpi::communicator{comm, boost::mpi::comm_take_ownership};
}

/** Returns a boost::mpi::communicator with a undirected Dist_graph
 * topology.
 */
inline boost::mpi::communicator
undirected_mpi_communicator(MPI_Comm parent_communicator,
                            const std::vector<int> &neighbor_ranks)
{
    return directed_mpi_communicator(parent_communicator, neighbor_ranks,
                                     neighbor_ranks);
}

} // namespace util
} // namespace repa
