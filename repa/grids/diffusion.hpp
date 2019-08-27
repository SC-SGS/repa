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
#include "glomethod.hpp"
#include "pargrid.hpp"
#include <array>
#include <mpi.h>
#include <unordered_map>
#include <vector>

namespace repa {
namespace grids {

/** Diffusively load-balanced grid.
 * Processes iteratively exchange boundary cells with neighbors.
 */
struct Diffusion : public GloMethod {
    Diffusion(const boost::mpi::communicator &comm,
              Vec3d box_size,
              double min_cell_size);
    ~Diffusion();

private:
    friend struct HybridGPDiff; // Needs access to "partition" vector

    void pre_init(bool firstcall) override;
    void post_init(bool firstcall) override;
    void init_new_foreign_cell(lidx localcell,
                               lgidx foreigncell,
                               rank owner) override;
    bool sub_repartition(CellMetric m, CellCellMetric ccm) override;

    //
    // Additional data
    //

    // Stores all local cells which ar at the border of the local area of
    // this process (Stores the local index)
    std::vector<int> borderCells;
    // Stores for each cell in "borderCells" the ranks of their neighbourhood
    // Key is the local cell ID and value a set of ranks
    std::map<lidx, std::vector<rank>> borderCellsNeighbors;

    // Neighborhood communicator
    MPI_Comm neighcomm;

    //
    // Additional methods
    //

    // Clears obsolete entries from "partition"
    void clear_unknown_cell_ownership();

    // Computes vector of vectors of cells which has to be send to neighbours
    std::vector<std::vector<int>>
    compute_send_list(std::vector<double> &&sendLoads,
                      const std::vector<double> &weights);

    // Struct for communication of neightbohood information
    struct NeighSend {
        int basecell;
        std::array<int, 26> neighranks;
    };

    // Send message with neighbourhood of received cells in "sendCells"
    std::vector<std::vector<NeighSend>>
    sendNeighbourhood(const std::vector<std::vector<int>> &toSend);

    // Update partition array
    void updateReceivedNeighbourhood(
        const std::vector<std::vector<NeighSend>> &neighbourhood);
};
} // namespace grids
} // namespace repa
