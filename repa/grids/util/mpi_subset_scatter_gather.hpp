
#pragma once

#include "pargrid.hpp"
#include <boost/mpi.hpp>

namespace repa {
namespace util {

/** MPI_Scatter "data" to the processes indicated by "neighbors" as well as
 * MPI_Gather from them and return this data.
 */
template <typename T>
std::vector<T>
mpi_subset_scatter_gather(const boost::mpi::communicator &comm,
                     const std::vector<grids::rank_type> &neighbors,
                     const std::vector<T> &data)
{
    std::vector<boost::mpi::request> sreq_cells(neighbors.size());
    std::vector<boost::mpi::request> rreq_cells(neighbors.size());

    for (grids::rank_index_type i = 0; i < neighbors.size(); ++i) {
        sreq_cells[i] = comm.isend(neighbors[i], 2, data[i]);
    }

    std::vector<T> gathered_data(neighbors.size());
    for (grids::rank_index_type i = 0; i < neighbors.size(); ++i) {
        rreq_cells[i] = comm.irecv(neighbors[i], 2, gathered_data[i]);
    }

    boost::mpi::wait_all(std::begin(rreq_cells), std::end(rreq_cells));
    boost::mpi::wait_all(std::begin(sreq_cells), std::end(sreq_cells));

    return gathered_data;
}


} // namespace util
} // namespace repa
