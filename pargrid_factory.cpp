
#include "pargrid_factory.hpp"
#include "grids/p4est.hpp"
#include "grids/cart.hpp"
#include "grids/graph.hpp"
#include "grids/kd_tree.hpp"
#include "grids/diffusion.hpp"
#include "grids/hybrid-gp-diff.hpp"
#include "grids/gridbased.hpp"

namespace generic_dd {
namespace grids {

namespace impl {
ParallelLCGrid *make_pargrid_impl(GridType gt) {
  switch (gt) {

  case GridType::P4EST:
#ifdef HAVE_P4EST
    return new P4estGrid();
#else
    throw std::invalid_argument("P4est not compiled in but requesting p4est grid.");
#endif

  case GridType::CART:
    return new CartGrid(); 

  case GridType::GRAPH:
#ifdef HAVE_METIS
    return new Graph();
#else
    throw std::invalid_argument("This ESPResSo has not been compiled with Metis support.");
#endif

  case GridType::DIFF:
    return new Diffusion();
    
  case GridType::KD_TREE:
#ifdef HAVE_KDPART
    return new KDTreeGrid();
#else
    throw std::invalid_argument("This ESPResSo has not been compiled with the kdpart-library.");
#endif

  case GridType::HYB_GP_DIFF:
#ifdef HAVE_METIS
    return new HybridGPDiff();
#else
    throw std::invalid_argument("This ESPResSo has not been compiled with Metis support.");
#endif

  case GridType::GRIDBASED:
#ifdef HAVE_TETRA
    return new GridBasedGrid();
#else
    throw std::invalid_argument("This ESPRsSo has not been compiled with Tetra support.");
#endif

  default:
    throw std::invalid_argument("Invalid grid type");
  }
}
}

std::unique_ptr<ParallelLCGrid>
make_pargrid(GridType gt) {
  return std::unique_ptr<ParallelLCGrid>(impl::make_pargrid_impl(gt));
}

}
}

