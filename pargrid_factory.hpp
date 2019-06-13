
#ifndef PARGRID_FACTORY_INCLUDED
#define  PARGRID_FACTORY_INCLUDED

#include <memory>
#include "pargrid.hpp"
#include "generic_dd_grid_types.hpp"

namespace generic_dd {
namespace grids {
/** Grid factory method.
 * To be called on every node.
 */
std::unique_ptr<ParallelLCGrid> make_pargrid(GridType gt);

}
}

#endif