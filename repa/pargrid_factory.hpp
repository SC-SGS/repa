
#pragma once

#include "grid_types.hpp"
#include "pargrid.hpp"
#include <memory>

namespace repa {
/** Grid factory method.
 * To be called on every node.
 */
std::unique_ptr<grids::ParallelLCGrid>
make_pargrid(GridType gt,
             const boost::mpi::communicator &comm,
             Vec3d box_size,
             double min_cell_size);

} // namespace repa
