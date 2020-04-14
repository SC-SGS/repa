/**
 * Copyright 2017-2019 The repa authors
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

#include "psdiffusion.hpp"

#include <algorithm>
#include <set>

#include "grids/util/mpi_cart.hpp"
#include "grids/util/mpi_graph.hpp"
#include "util/vec_arith.hpp"

#ifndef NDEBUG
#define PSDIFFUSION_DEBUG
#endif

namespace repa {
namespace grids {

/** Returns the new Coords of a Rank after it beeing maped to opposite site.
 * This is a necessary step for the metric because of the periodic edge.
 */
static Vec3i map_coords_to_opposite_side(const Vec3i &c0,
                                         const Vec3i &c2,
                                         const Vec3i &gsize)
{
    using namespace util::vector_arithmetic;
    Vec3i te, ts;
    for (int i = 0; i < te.size(); i++) {
        te[i] = c2[i] - c0[i] > 1;
        ts[i] = c0[i] - c2[i] > 1;
    }
    return c2 + ts * gsize - te * gsize;
}

/*
 * Initialization
 */
PSDiffusion::PSDiffusion(const boost::mpi::communicator &comm,
                         Vec3d box_size,
                         double min_cell_size)
    : Diffusion(comm, box_size, min_cell_size)
{
}

PSDiffusion::~PSDiffusion()
{
}

void PSDiffusion::post_init(bool firstcall)
{
    Diffusion::post_init(firstcall);
    comm_dims = util::mpi_cart_get_dims(comm_cart);

#ifdef PSDIFFUSION_DEBUG
    // Copy initial neighborhood, so that the neighborhood can be checked for
    // consistency in later iterations.
    if (firstcall) {
        initial_neighborhood
            = std::vector<rank_type>(neighbors.begin(), neighbors.end());
    }

    // Neighborhood consistency check
    std::set<rank_type> a(initial_neighborhood.begin(),
                          initial_neighborhood.end());
    std::set<rank_type> b(neighbors.begin(), neighbors.end());

    assert(std::includes(a.begin(), a.end(), b.begin(), b.end()));
#endif
}

bool PSDiffusion::accept_transfer(local_cell_index_type cidx,
                                  rank_type neighrank) const
{
    return coords_based_allow_sending(cidx, neighrank);
}

bool PSDiffusion::coords_based_allow_sending(local_cell_index_type c,
                                             rank_type neighrank) const
{
    Vec3i c0, cn;
    MPI_Cart_coords(comm_cart, rank_of_cell(cells[c]), 3, c0.data());
    MPI_Cart_coords(comm_cart, neighrank, 3, cn.data());
    for (const global_cell_index_type &d :
         gbox.full_shell_neigh_without_center(cells[c])) {
        rank_type rank_d = rank_of_cell(d);
        if (rank_d == rank_of_cell(cells[c]) || rank_d == neighrank)
            continue;
        Vec3i c2;
        MPI_Cart_coords(comm_cart, rank_d, 3, c2.data());

        cn = map_coords_to_opposite_side(c0, cn, comm_dims);
        c2 = map_coords_to_opposite_side(c0, c2, comm_dims);
        if (std::abs(cn[0] - c2[0]) > 2 || std::abs(cn[1] - c2[1]) > 2
            || std::abs(cn[2] - c2[2]) > 2)
            return false;
    }
    return true;
}
} // namespace grids
} // namespace repa
