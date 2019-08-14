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
#include "grids/graph.hpp"
#include "grids/gridbased.hpp"
#include "grids/hybrid-gp-diff.hpp"
#include "grids/kd_tree.hpp"
#include "grids/p4est.hpp"

namespace repa {
namespace grids {

namespace {
ParallelLCGrid *make_pargrid_impl(GridType gt,
                                  const boost::mpi::communicator &comm,
                                  Vec3d box_size,
                                  double min_cell_size,
                                  ExtraParams ep)
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
        return new GridBasedGrid(comm, box_size, min_cell_size, ep);
        //#else
        //    throw std::invalid_argument("This ESPRsSo has not been compiled
        //    with Tetra support.");
        //#endif

    default:
        throw std::invalid_argument("Invalid grid type");
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
