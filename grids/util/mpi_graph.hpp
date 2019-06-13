
#pragma once

#include <mpi.h>
#include "common_types.hpp"

namespace repa {
namespace util {

/** Returns the number of neighbors in an undirected Dist_graph communicator.
 * This function simply assumes that the graph topology is undirected and
 * does not verify it!
 */
inline int mpi_undirected_neighbor_count(MPI_Comm neighcomm) {
  int indegree = 0, outdegree = 0, weighted = 0;
  MPI_Dist_graph_neighbors_count(neighcomm, &indegree, &outdegree, &weighted);
  return indegree;
}

}
}
