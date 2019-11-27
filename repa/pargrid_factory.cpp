/**
 * Copyright 2017-2019 Steffen Hirschmann
 *
 * This file is part of Repa.
 *
 * Repa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Repa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Repa.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pargrid_factory.hpp"
#include "grids/cart.hpp"
#include "grids/diffusion.hpp"
#include "grids/psdiffusion.hpp"
#ifdef HAVE_PARMETIS
#include "grids/graph.hpp"
#include "grids/hybrid-gp-diff.hpp"
#endif
#ifdef HAVE_CGAL
#include "grids/gridbased.hpp"
#endif
#ifdef HAVE_KDPART
#include "grids/kd_tree.hpp"
#endif
#ifdef HAVE_P4EST
#include "grids/p4est.hpp"
#endif

namespace repa {
namespace grids {

namespace {
ParallelLCGrid *make_pargrid_impl(GridType gt,
                                  const boost::mpi::communicator &comm,
                                  Vec3d box_size,
                                  double min_cell_size,
                                  ExtraParams ep)
{
    ParallelLCGrid *r = nullptr;
    switch (gt) {
    case GridType::P4EST:
#ifdef HAVE_P4EST
        r = new P4estGrid(comm, box_size, min_cell_size);
#else
        throw std::invalid_argument("Librepa not compiled with p4est support.");
#endif
        break;

    case GridType::CART:
        r = new CartGrid(comm, box_size, min_cell_size);
        break;

    case GridType::GRAPH:
#ifdef HAVE_PARMETIS
        r = new Graph(comm, box_size, min_cell_size);
#else
        throw std::invalid_argument(
            "Librepa not compiled with Parmetis support.");
#endif
        break;

    case GridType::DIFF:
        r = new Diffusion(comm, box_size, min_cell_size);
        break;

    case GridType::PSDIFF:
        r = new PSDiffusion(comm, box_size, min_cell_size);
        break;

    case GridType::KD_TREE:
#ifdef HAVE_KDPART
        r = new KDTreeGrid(comm, box_size, min_cell_size);
#else
        throw std::invalid_argument(
            "Librepa not compiled with kdpart support.");
#endif
        break;

    case GridType::HYB_GP_DIFF:
#ifdef HAVE_PARMETIS
        r = new HybridGPDiff(comm, box_size, min_cell_size);
#else
        throw std::invalid_argument(
            "Librepa not compiled with Parmetis support.");
#endif
        break;

    case GridType::GRIDBASED:
#ifdef HAVE_CGAL
        r = new GridBasedGrid(comm, box_size, min_cell_size, ep);
#else
        throw std::invalid_argument("Librepa not compiled with CGAL support");
#endif
        break;

    default:
        throw UnknownGridTypeError(std::to_string(static_cast<int>(gt)));
        break;
    }

    if (r) {
        r->after_construction();
        return r;
    }
    else {
        throw std::runtime_error("Did not construct a grid?");
    }
}
} // namespace
} // namespace grids

std::unique_ptr<grids::ParallelLCGrid>
make_pargrid(GridType gt,
             const boost::mpi::communicator &comm,
             Vec3d box_size,
             double min_cell_size,
             ExtraParams ep)
{
    return std::unique_ptr<grids::ParallelLCGrid>(
        grids::make_pargrid_impl(gt, comm, box_size, min_cell_size, ep));
}

} // namespace repa
