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
#include "pargrid.hpp" // UNKNOWN_RANK
#include <boost/mpi/communicator.hpp>
#include <mpi.h>

namespace repa {
namespace util {

/** Returns true if "comm" has DIST_GRAPH topology.
 */
inline bool has_dist_graph_topology(MPI_Comm comm)
{
    int status;
    MPI_Topo_test(comm, &status);
    return status == MPI_DIST_GRAPH;
}

/** Returns the number of neighbors in an undirected Dist_graph communicator.
 * This function simply assumes that the graph topology is undirected and
 * does not verify it!
 */
inline int mpi_undirected_neighbor_count(MPI_Comm neighcomm)
{
    assert(has_dist_graph_topology(neighcomm));

    int indegree = 0, outdegree = 0, weighted = 0;
    MPI_Dist_graph_neighbors_count(neighcomm, &indegree, &outdegree, &weighted);
    assert(indegree == outdegree);
    return indegree;
}

inline std::vector<int> mpi_undirected_neighbors(MPI_Comm neighcomm)
{
    assert(has_dist_graph_topology(neighcomm));

    auto nneigh = util::mpi_undirected_neighbor_count(neighcomm);

    std::vector<int> srcneigh(nneigh, -1),
        dstneigh(nneigh, -1), dummy(nneigh);
    MPI_Dist_graph_neighbors(neighcomm, nneigh, srcneigh.data(), dummy.data(),
                             nneigh, dstneigh.data(), dummy.data());

    assert(srcneigh == dstneigh); // Undirected sanity check
    assert(std::find(srcneigh.begin(), srcneigh.end(), -1)
           == srcneigh.end());
    return srcneigh;
}

inline std::pair<std::vector<int>, std::vector<int>>
mpi_directed_neighbors(MPI_Comm neighcomm)
{
    assert(has_dist_graph_topology(neighcomm));

    int indegree = 0, outdegree = 0, weighted = 0;
    MPI_Dist_graph_neighbors_count(neighcomm, &indegree, &outdegree, &weighted);

    std::vector<int> srcneigh(indegree, -1), dummy1(indegree),
        dstneigh(outdegree, -1), dummy2(outdegree);
    MPI_Dist_graph_neighbors(neighcomm, indegree, srcneigh.data(),
                             dummy1.data(), outdegree, dstneigh.data(),
                             dummy2.data());

    assert(std::find(srcneigh.begin(), srcneigh.end(), -1)
           == srcneigh.end());
    assert(std::find(dstneigh.begin(), dstneigh.end(), -1)
           == dstneigh.end());
    return std::make_pair(srcneigh, dstneigh);
}

/** Returns a boost::mpi::communicator with a undirected Dist_graph
 * topology.
 */
inline boost::mpi::communicator
directed_graph_communicator(MPI_Comm parent_communicator,
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
undirected_graph_communicator(MPI_Comm parent_communicator,
                              const std::vector<int> &neighbor_ranks)
{
    return directed_graph_communicator(parent_communicator, neighbor_ranks,
                                       neighbor_ranks);
}

} // namespace util
} // namespace repa
