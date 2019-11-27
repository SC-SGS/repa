/**
 * Copyright 2017-2019 Steffen Hirschmann
 * Copyright 2017-2018 Maximilian Wildbrett
 * Copyright 2019      Simon Hauser
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

#include "grids/util/mpi_graph.hpp"
#include "util/ensure.hpp"

#ifndef NDEBUG
#define PSDIFFUSION_DEBUG
#endif

namespace repa {
namespace grids {

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

    std::vector<rank_type> send_ranks(neighbors.begin(), neighbors.end());
    send_ranks.insert(send_ranks.begin(), comm_cart.rank());

    neighborhood_ranks.resize(neighbors.size() * send_ranks.size());
    MPI_Neighbor_allgather(send_ranks.data(), send_ranks.size(), MPI_INT,
                           neighborhood_ranks.data(), send_ranks.size(),
                           MPI_INT, neighcomm);
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
        ENSURE(nadditional_comm < 27);
#endif

        if (profit > 0)
            plist.emplace_back(27 - nadditional_comm, profit, borderCells[i]);
    }

    PerNeighbor<GlobalCellIndices> to_send(send_loads.size());

    // Use a maxheap: Always draw the maximum element
    // (1. least new border cells, 2. most profit)
    // and find a process that can take this cell.
    std::make_heap(std::begin(plist), std::end(plist));
#ifdef PSDIFFUSION_DEBUG
    int droppedCells = 0;
    int maxCells = 0;
#endif
    while (!plist.empty()) {
        std::pop_heap(std::begin(plist), std::end(plist));
        local_cell_index_type cidx = std::get<2>(plist.back());
        plist.pop_back();
#ifdef PSDIFFUSION_DEBUG
        maxCells++;
#endif

        for (auto neighrank : borderCellsNeighbors[cidx]) {
            auto neighidx
                = std::distance(std::begin(neighbors),
                                std::find(std::begin(neighbors),
                                          std::end(neighbors), neighrank));

            bool b1 = coords_based_allow_sending(cidx, neighrank);
#ifdef PSDIFFUSION_DEBUG
            bool b2 = rank_based_allow_sending(cidx, neighrank);
            if (!b1)
                droppedCells++;
            if (!b2)
                std::cout << "Doing work" << std::endl;

                /* if (b1 != b2) */
                /*     std::cout << "B1: " << b1 << "; B2: " << b2 << std::endl;
                 */
#endif
            if (!b1)
                continue;

            if (weights[cidx] <= send_loads[neighidx]) {
                to_send[neighidx].push_back(cells[cidx]);
                send_loads[neighidx] -= weights[cidx];
                // This cell is done. Continue with the next.
                break;
            }
        }
    }
#ifdef PSDIFFUSION_DEBUG
    std::cout << "Dropped: " << droppedCells << "/" << maxCells << " Cells"
              << std::endl;
#endif

    return to_send;
}

bool PSDiffusion::rank_based_allow_sending(local_cell_index_type c,
                                           rank_type neighrank)
{
    for (global_cell_index_type d1 :
         gbox.full_shell_neigh_without_center(cells[c])) {
        if (rank_of_cell(d1) != neighrank)
            continue;
        for (global_cell_index_type d2 :
             gbox.full_shell_neigh_without_center(cells[c])) {
            rank_type r1 = rank_of_cell(d1);
            rank_type r2 = rank_of_cell(d2);
            if (r1 == r2 || r2 == rank_of_cell(cells[c]))
                continue;
#ifdef PSDIFFUSION_DEBUG
            ENSURE(neighborhood_ranks.size() != 0);
#endif

            bool found = false;
            for (int y = 0; y < neighbors.size(); y++) {
                for (int x = 0; x < neighbors.size() + 1; x++) {
                    int index = x + (neighbors.size() + 1) * y;
                    if (x == 0) {
                        if (neighborhood_ranks[index] != r2)
                            break;
                    }
                    else {
                        if (neighborhood_ranks[index] == r1)
                            found = true;
                    }
                }
            }
            if (!found)
                return false;
        }
    }
    return true;
}

bool PSDiffusion::coords_based_allow_sending(local_cell_index_type c,
                                             rank_type neighrank)
{
    for (global_cell_index_type d1 :
         gbox.full_shell_neigh_without_center(cells[c])) {
        if (rank_of_cell(d1) != neighrank)
            continue;

        Vec3i c1;
        MPI_Cart_coords(comm_cart, rank_of_cell(d1), 3, c1.data());
        for (global_cell_index_type d2 :
             gbox.full_shell_neigh_without_center(cells[c])) {

            Vec3i c2;
            MPI_Cart_coords(comm_cart, rank_of_cell(d2), 3, c2.data());

            if (std::abs(c1[0] - c2[0]) > 2 || std::abs(c1[1] - c2[1]) > 2
                || std::abs(c1[2] - c2[2]) > 2)
                return false;
        }
    }
    return true;
}

} // namespace grids
} // namespace repa
