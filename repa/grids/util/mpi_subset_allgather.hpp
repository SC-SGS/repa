
#pragma once

#include "pargrid.hpp"
#include <boost/mpi.hpp>

namespace repa {
namespace util {

/** MPI_Allgather to and from a subset of processes
 *
 */
template <typename T, typename U = T>
std::vector<T>
mpi_subset_allgather(const boost::mpi::communicator &comm,
                     const std::vector<grids::rank_type> &neighbors,
                     const U &data)
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

namespace __impl {

/**
 * Type trait that has a boolean member "value" equal to true if the
 * template parameter is a std::pair. Otherwise, the member is false.
 */
template <typename T>
struct is_pair {
    static const bool value = false;
};

template <typename T, typename U>
struct is_pair<std::pair<T, U>> {
    static const bool value = true;
};

template <typename T>
const bool is_pair_v = is_pair<T>::value;
}

/** MPI_Allgather to and from a subset of processes
 * 
 * Overload specifically to send std::pairs that hold references.
 * Pass references and get out non-references.
 */
template <typename T1,
          typename T2,
          // Only enabled if "T1" is not a std::pair
          typename = typename std::enable_if<!__impl::is_pair_v<T1>>::type>
std::vector<std::pair<T1, T2>>
mpi_subset_allgather(const boost::mpi::communicator &comm,
                     const std::vector<grids::rank_type> &neighbors,
                     const std::pair<const T1 &, const T2 &> &data)
{
    return mpi_subset_allgather<std::pair<T1, T2>,
                                std::pair<const T1 &, const T2 &>>(
        comm, neighbors, data);
}

} // namespace util
} // namespace repa
