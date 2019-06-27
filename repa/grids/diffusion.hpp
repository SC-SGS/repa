/**
 * Copyright 2017-2019 Steffen Hirschmann
 * Copyright 2017-2018 Maximilian Wildbrett
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

#pragma once

#include "globox.hpp"
#include "pargrid.hpp"
#include <array>
#include <mpi.h>
#include <unordered_map>
#include <vector>

namespace repa {
namespace grids {

struct Diffusion : public ParallelLCGrid {
    Diffusion(const boost::mpi::communicator &comm,
              Vec3d box_size,
              double min_cell_size);
    ~Diffusion();
    lidx n_local_cells() override;
    gidx n_ghost_cells() override;
    nidx n_neighbors() override;
    rank neighbor_rank(nidx i) override;
    Vec3d cell_size() override;
    Vec3i grid_size() override;
    lgidx cell_neighbor_index(lidx cellidx, int neigh) override;
    std::vector<GhostExchangeDesc> get_boundary_info() override;
    lidx position_to_cell_index(double pos[3]) override;
    rank position_to_rank(double pos[3]) override;
    nidx position_to_neighidx(double pos[3]) override;
    bool repartition(CellMetric m,
                     CellCellMetric ccm,
                     Thunk exchange_start_callback) override;

    int global_hash(lgidx cellidx) override;

    struct NeighSend {
        int basecell;
        std::array<int, 26> neighranks;
    };

private:
    friend struct HybridGPDiff; // Needs access to "partition" vector
    // Number of local cells of this process(partition)
    lidx localCells;
    // Number of ghost cells of this process(partition)
    int ghostCells;
    // Vector with all neighbour ranks
    std::vector<rank> neighbors;
    // Stores the "GhostExchangeDesc" elements of the process(partition)
    std::vector<GhostExchangeDesc> exchangeVector;

    globox::GlobalBox<int, int> gbox;

    // Stores the local cells and the ghost cells of this process(partition)
    std::vector<int> cells;
    // Stores all local cells which ar at the border of the local area of
    // this process (Stores the local index)
    std::vector<int> borderCells;
    // Stores for each cell in "borderCells" the ranks of their neighbourhood
    // Key is the local cell ID and value a set of ranks
    std::map<lidx, std::vector<rank>> borderCellsNeighbors;
    // Stores the received partition array
    // Can be dircetly used in Metis algorithm
    std::vector<rank> partition;
    // Stores the mapping of a global linearized index to a
    // local linearized index or an ghost linearized index.
    // Global index is key and local/ghost index is value.
    std::unordered_map<int, int> global_to_local;

    // After partition rebuild the local part of the structure
    void reinit(bool init = false);

    // Clears obsolete entries from "partition"
    void clear_unknown_cell_ownership();

    // Computes vector of vectors of cells which has to be send to neighbours
    std::vector<std::vector<int>>
    compute_send_list(std::vector<double> &&sendLoads,
                      const std::vector<double> &weights);

    MPI_Comm neighcomm;

    // Send message with neighbourhood of received cells in "sendCells"
    std::vector<std::vector<NeighSend>>
    sendNeighbourhood(const std::vector<std::vector<int>> &toSend);
    // Update partition array
    void updateReceivedNeighbourhood(
        const std::vector<std::vector<NeighSend>> &neighbourhood);
};
} // namespace grids
} // namespace repa
