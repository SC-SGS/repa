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

static Vec3i
fix_periodic_edge(const Vec3i &c0, const Vec3i &c2, const Vec3i &gsize)
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
    if (firstcall) {
        initial_neighborhood
            = std::vector<rank_type>(neighbors.begin(), neighbors.end());

        // Allgather over neighbors.size
        std::vector<int> neighbors_sizes(neighbors.size());
        int nsize = neighbors.size();
        MPI_Neighbor_allgather(&nsize, 1, MPI_INT, neighbors_sizes.data(), 1,
                               MPI_INT, neighcomm);
        // Assert neighbors.size equal
        assert(std::equal(neighbors_sizes.begin() + 1, neighbors_sizes.end(),
                          neighbors_sizes.begin()));

        neighborhood_ranks.resize(neighbors.size() * neighbors.size());
        MPI_Neighbor_allgather(neighbors.data(), neighbors.size(), MPI_INT,
                               neighborhood_ranks.data(), neighbors.size(),
                               MPI_INT, neighcomm);

        for (int i = 0; i < neighbors.size(); i++)
            nr_mappings[neighbors[i]]
                = &neighborhood_ranks.data()[0] + (neighbors.size() * i);
    }

    // Neighborhood consistence check
    std::set<rank_type> a(initial_neighborhood.begin(),
                          initial_neighborhood.end());
    std::set<rank_type> b(neighbors.begin(), neighbors.end());

    assert(std::includes(a.begin(), a.end(), b.begin(), b.end()));
#endif
}

bool PSDiffusion::accept_transfer(local_cell_index_type cidx,
                                  rank_type neighrank)
{
    const bool b1 = coords_based_allow_sending(cidx, neighrank);
#ifdef PSDIFFUSION_DEBUG
    assert(b1 == rank_based_allow_sending(cidx, neighrank));
#endif
    return b1;
}

bool PSDiffusion::coords_based_allow_sending(local_cell_index_type c,
                                             rank_type neighrank)
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

        cn = fix_periodic_edge(c0, cn, comm_dims);
        c2 = fix_periodic_edge(c0, c2, comm_dims);
        if (std::abs(cn[0] - c2[0]) > 2 || std::abs(cn[1] - c2[1]) > 2
            || std::abs(cn[2] - c2[2]) > 2)
            return false;
    }
    return true;
}

#ifdef PSDIFFUSION_DEBUG
bool PSDiffusion::rank_based_allow_sending(local_cell_index_type c,
                                           rank_type neighrank)
{
    for (const global_cell_index_type &d1 :
         gbox.full_shell_neigh_without_center(cells[c])) {
        if (rank_of_cell(d1) != neighrank)
            continue;
        rank_type r1 = rank_of_cell(d1);
        auto npr_begin = nr_mappings[r1];
        auto npr_end = npr_begin + initial_neighborhood.size();
        for (const global_cell_index_type &d2 :
             gbox.full_shell_neigh(cells[c])) {
            rank_type r2 = rank_of_cell(d2);
            if (r1 == r2 || r2 == rank_of_cell(cells[c]))
                continue;
            if (std::find(npr_begin, npr_end, r2) == npr_end)
                return false;
        }
    }
    return true;
}
#endif

} // namespace grids
} // namespace repa
