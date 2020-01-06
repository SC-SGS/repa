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

#include <set>

#include "grids/util/mpi_graph.hpp"
#include "util/vec_arith.hpp"

#ifndef NDEBUG
#define PSDIFFUSION_DEBUG
#endif

#ifdef PSDIFFUSION_DEBUG
/*
 * @returns True if structure is fine, false if not
 */
static bool check_structure(std::vector<repa::grids::rank_type> initial,
                            std::vector<repa::grids::rank_type> current)
{
    std::set<repa::grids::rank_type> L1(initial.begin(), initial.end());
    std::set<repa::grids::rank_type> L2(current.begin(), current.end());
    auto it1 = L1.begin();
    auto it2 = L2.begin();

    while (it1 != L1.end() && it2 != L2.end()) {
        if (*it1 < *it2) {
            // Losing Neighbors is allowed
            *it1++;
        }
        else if (*it2 < *it1) {
            // Gaining Neighbors is not allowed
            return false;
        }
        else {
            it1++;
            it2++;
        }
    }

    // Gained at least 1 Neighbor so true
    if (it2 != L2.end())
        return false;

    return true;
}
#endif

namespace repa {
namespace grids {

std::vector<double> PSDiffusion::compute_send_volume(double load)
{
#ifdef PSDIFFUSION_DEBUG
    assert(check_structure(initial_neighborhood, neighbors));
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
    std::vector<rank_type> send_ranks(neighbors.begin(), neighbors.end());
    send_ranks.insert(send_ranks.begin(), comm_cart.rank());

    neighborhood_ranks.resize(neighbors.size() * send_ranks.size());
    MPI_Neighbor_allgather(send_ranks.data(), send_ranks.size(), MPI_INT,
                           neighborhood_ranks.data(), send_ranks.size(),
                           MPI_INT, neighcomm);

    if (firstcall) {
        initial_neighborhood
            = std::vector<rank_type>(neighbors.begin(), neighbors.end());
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

            if (b1 != b2)
                std::cout << "B1: " << b1 << "; B2: " << b2 << std::endl;
            /* assert(b1 == b2); */
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

    return c2 + ts * 4 - te * 4;
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
        auto npr = neighbar_procs(r1);
        for (const global_cell_index_type &d2 :
             gbox.full_shell_neigh(cells[c])) {
            rank_type r2 = rank_of_cell(d2);
            if (r1 == r2)
                continue;
            if (std::find(npr.begin(), npr.end(), r2) == npr.end())
                return false;
        }
    }
    return true;
}

std::vector<rank_type> PSDiffusion::neighbar_procs(rank_type r)
{
    int startingIndex = -1;
    int endingIndex = -1;
    for (int y = 0; y < neighbors.size(); y++) {
        int index = (neighbors.size() + 1) * y;
        if (neighborhood_ranks[index] == r) {
            startingIndex = index + 1;
            endingIndex = startingIndex + neighbors.size();
        }
    }

    assert(startingIndex != -1 && endingIndex != -1);

    auto start = neighborhood_ranks.begin() + startingIndex;
    auto end = neighborhood_ranks.begin() + endingIndex;

    auto result = std::vector<rank_type>(start, end);

    assert(result.size() == neighbors.size());

    return std::vector<rank_type>(start, end);
}
#endif

} // namespace grids
} // namespace repa
