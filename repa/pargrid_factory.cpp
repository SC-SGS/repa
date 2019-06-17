
#include "pargrid_factory.hpp"
#include "grids/cart.hpp"
#include "grids/diffusion.hpp"
#include "grids/graph.hpp"
#include "grids/gridbased.hpp"
#include "grids/hybrid-gp-diff.hpp"
#include "grids/kd_tree.hpp"
#include "grids/p4est.hpp"

namespace repa {
namespace grids {

namespace impl {
ParallelLCGrid *make_pargrid_impl(GridType gt,
                                  const boost::mpi::communicator &comm,
                                  Vec3d box_size,
                                  double min_cell_size)
{
    switch (gt) {

    case GridType::P4EST:
        //#ifdef HAVE_P4EST
        return new P4estGrid(comm, box_size, min_cell_size);
        //#else
        //    throw std::invalid_argument("P4est not compiled in but requesting
        //    p4est grid.");
        //#endif

    case GridType::CART:
        return new CartGrid(comm, box_size, min_cell_size);

    case GridType::GRAPH:
        //#ifdef HAVE_METIS
        return new Graph(comm, box_size, min_cell_size);
        //#else
        //    throw std::invalid_argument("This ESPResSo has not been compiled
        //    with Metis support.");
        //#endif

    case GridType::DIFF:
        return new Diffusion(comm, box_size, min_cell_size);

    case GridType::KD_TREE:
        //#ifdef HAVE_KDPART
        return new KDTreeGrid(comm, box_size, min_cell_size);
        //#else
        //    throw std::invalid_argument("This ESPResSo has not been compiled
        //    with the kdpart-library.");
        //#endif

    case GridType::HYB_GP_DIFF:
        //#ifdef HAVE_METIS
        return new HybridGPDiff(comm, box_size, min_cell_size);
        //#else
        //    throw std::invalid_argument("This ESPResSo has not been compiled
        //    with Metis support.");
        //#endif

    case GridType::GRIDBASED:
        //#ifdef HAVE_TETRA
        return new GridBasedGrid(comm, box_size, min_cell_size);
        //#else
        //    throw std::invalid_argument("This ESPRsSo has not been compiled
        //    with Tetra support.");
        //#endif

    default:
        throw std::invalid_argument("Invalid grid type");
    }
}
} // namespace impl

std::unique_ptr<ParallelLCGrid>
make_pargrid(GridType gt,
             const boost::mpi::communicator &comm,
             Vec3d box_size,
             double min_cell_size)
{
    return std::unique_ptr<ParallelLCGrid>(
        impl::make_pargrid_impl(gt, comm, box_size, min_cell_size));
}

} // namespace grids
} // namespace repa
