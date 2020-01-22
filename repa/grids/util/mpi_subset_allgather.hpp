
#pragma once

#include "pargrid.hpp"
#include <boost/mpi.hpp>

namespace repa {
namespace util {

template <typename T>
std::vector<T>
mpi_subset_allgather(const boost::mpi::communicator &comm,
                     const std::vector<grids::rank_type> &neighbors,
                     const T &data)
{
    std::vector<boost::mpi::request> sreq_cells(neighbors.size());
    std::vector<boost::mpi::request> rreq_cells(neighbors.size());

    // Extra loop as all ranks need to be added before sending
    for (grids::rank_index_type i = 0; i < neighbors.size(); ++i) {
        sreq_cells[i] = comm.isend(neighbors[i], 2, data);
    }

    // All send volumes from all processes
    std::vector<T> all_data(neighbors.size());
    for (grids::rank_index_type i = 0; i < neighbors.size(); ++i) {
        rreq_cells[i] = comm.irecv(neighbors[i], 2, all_data[i]);
    }

    boost::mpi::wait_all(std::begin(rreq_cells), std::end(rreq_cells));
    boost::mpi::wait_all(std::begin(sreq_cells), std::end(sreq_cells));

    return all_data;
}

} // namespace util
} // namespace repa
