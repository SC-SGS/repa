
#include "pargrid.hpp"

#include "grids/util/mpi_cart.hpp"

namespace repa {
namespace grids {

static boost::mpi::communicator
make_cartesian_communicator(const boost::mpi::communicator &comm)
{
    Vec3i dims = {{0, 0, 0}}, periods = {{1, 1, 1}};
    MPI_Dims_create(comm.size(), 3, dims.data());

    MPI_Comm _comm_cart;
    MPI_Cart_create(comm, 3, dims.data(), periods.data(), true, &_comm_cart);

    return boost::mpi::communicator{_comm_cart,
                                    boost::mpi::comm_take_ownership};
}

ParallelLCGrid::ParallelLCGrid(const boost::mpi::communicator &comm,
                               Vec3d box_size,
                               double min_cell_size)
    : comm(comm, boost::mpi::comm_duplicate),
      comm_cart(make_cartesian_communicator(comm)),
      box_l(box_size),
      node_grid(util::mpi_cart_get_dims(comm_cart)),
      node_pos(util::mpi_cart_get_coords(comm_cart)),
      max_range(min_cell_size)
{
}

} // namespace grids
} // namespace repa