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

} // namespace util
} // namespace repa
