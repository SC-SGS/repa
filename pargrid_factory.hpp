
#pragma once

#include "generic_dd_grid_types.hpp"
#include "pargrid.hpp"
#include <memory>

namespace repa {
namespace grids {
/** Grid factory method.
 * To be called on every node.
 */
std::unique_ptr<ParallelLCGrid>
make_pargrid(GridType gt,
             const boost::mpi::communicator &comm,
             Vec3d box_size,
             double min_cell_size);

} // namespace grids
} // namespace repa
