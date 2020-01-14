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

std::vector<double> PSDiffusion::compute_send_volume(double load)
{
#ifdef PSDIFFUSION_DEBUG
    std::set<rank_type> a(initial_neighborhood.begin(),
                          initial_neighborhood.end());
    std::set<rank_type> b(neighbors.begin(), neighbors.end());

    assert(std::includes(a.begin(), a.end(), b.begin(), b.end()));
#endif

    return Diffusion::compute_send_volume(load);
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
#endif
}

/*
 * Computes a vector of vectors. The inner vectors contain a rank of the
 * process where the cells shall send and the cellids of this cells.
 */
Diffusion::PerNeighbor<Diffusion::GlobalCellIndices>
PSDiffusion::compute_send_list(std::vector<double> &&send_loads,
                               const std::vector<double> &weights)
{
    std::vector<std::tuple<int, double, local_cell_index_type>> plist;
    for (size_t i = 0; i < borderCells.size(); i++) {
        // Profit when sending this cell away
        double profit = weights[borderCells[i]];

        // Additional cell communication induced if this cell is sent away
        int nadditional_comm = 0;
        for (global_cell_index_type neighCell :
             gbox.full_shell_neigh_without_center(cells[borderCells[i]])) {
            if (partition[neighCell] == comm_cart.rank()
                && std::find(std::begin(borderCells), std::end(borderCells),
                             global_to_local[neighCell])
                       != std::end(borderCells)) {
                nadditional_comm++;
            }
        }
#ifdef DIFFUSION_DEBUG
        assert(nadditional_comm < 27);
#endif

        if (profit > 0)
            plist.emplace_back(27 - nadditional_comm, profit, borderCells[i]);
    }

    PerNeighbor<GlobalCellIndices> to_send(send_loads.size());

    // Use a maxheap: Always draw the maximum element
    // (1. least new border cells, 2. most profit)
    // and find a process that can take this cell.
    std::make_heap(std::begin(plist), std::end(plist));
    while (!plist.empty()) {
        std::pop_heap(std::begin(plist), std::end(plist));
        local_cell_index_type cidx = std::get<2>(plist.back());
        plist.pop_back();

        for (auto neighrank : borderCellsNeighbors[cidx]) {
            bool b1 = coords_based_allow_sending(cidx, neighrank);
#ifdef PSDIFFUSION_DEBUG
            bool b2 = rank_based_allow_sending(cidx, neighrank);
            assert(b1 == b2);
#endif
            if (!b1)
                continue;

            auto neighidx
                = std::distance(std::begin(neighbors),
                                std::find(std::begin(neighbors),
                                          std::end(neighbors), neighrank));

            if (weights[cidx] <= send_loads[neighidx]) {
                to_send[neighidx].push_back(cells[cidx]);
                send_loads[neighidx] -= weights[cidx];
                // This cell is done. Continue with the next.
                break;
            }
        }
    }

    return to_send;
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

        cn = fix_periodic_edge(c0, cn);
        c2 = fix_periodic_edge(c0, c2);
        if (std::abs(cn[0] - c2[0]) > 2 || std::abs(cn[1] - c2[1]) > 2
            || std::abs(cn[2] - c2[2]) > 2)
            return false;
    }
    return true;
}

Vec3i PSDiffusion::fix_periodic_edge(const Vec3i &c0, const Vec3i &c2)
{
    using namespace util::vector_arithmetic;
    Vec3i te, ts;
    for (int i = 0; i < te.size(); i++) {
        te[i] = c2[i] - c0[i] > 1;
        ts[i] = c0[i] - c2[i] > 1;
    }

    Vec3i gsize = util::mpi_cart_get_dims(comm_cart);

    return c2 + ts * gsize - te * gsize;
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
