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

#pragma once

#include "diff_variants.hpp"
#include "globox.hpp"
#include "glomethod.hpp"
#include "pargrid.hpp"
#include <array>
#include <set>
#include <vector>

namespace repa {
namespace grids {

namespace __diff_impl {
// Struct for communication of neightbohood information
struct CellNeighborhood {
    global_cell_index_type basecell;
    std::array<rank_type, 26> neighranks;
};

using CellNeighborhoodPerCell = std::vector<CellNeighborhood>;
} // namespace __diff_impl

/** Diffusively load-balanced grid.
 * Processes iteratively exchange boundary cells with neighbors.
 */
struct Diffusion : public GloMethod {
    Diffusion(const boost::mpi::communicator &comm,
              Vec3d box_size,
              double min_cell_size);
    ~Diffusion();

    virtual std::set<std::string> get_supported_variants();

protected:
    virtual void post_init(bool firstcall) override;

    virtual rank_type rank_of_cell(global_cell_index_type idx) const override
    {
        assert(idx >= 0 && idx < gbox.ncells());
        return partition[idx];
    }

    //
    // Additional data
    //

    // Stores all local cells which ar at the border of the local area of
    // this process (Stores the local index)
    std::vector<local_cell_index_type> borderCells;

    // Stores for each cell in "borderCells" the ranks of their neighborhood
    // Key is the local cell ID and value a set of ranks
    std::map<local_cell_index_type, std::vector<rank_type>>
        borderCellsNeighbors;

    // Stores the global partitioning. One rank per cell. Index via global
    // index.
    std::vector<rank_type> partition;

    // Neighborhood communicator
    boost::mpi::communicator neighcomm;

    /** Type "Per_Neighbor": an element designated for communication with
     * a neighboring process (one per rank_index_type).
     */
    template <typename T>
    using PerNeighbor = std::vector<T>;

    /** Communication volume (list of cells)
     */
    using GlobalCellIndices = std::vector<global_cell_index_type>;

    virtual void command(std::string s) override;

    virtual bool accept_transfer(local_cell_index_type cell,
                                 rank_type rank) const
    {
        return true;
    };

    std::unique_ptr<diff_variants::FlowCalculator> flow_calc;

private:
    friend struct HybridGPDiff; // Needs access to "partition" vector

    double profit_percentage_pass_through = 1.0;

    void pre_init(bool firstcall) override;
    void init_new_foreign_cell(local_cell_index_type localcell,
                               global_cell_index_type foreigncell,
                               rank_type owner) override;

    bool sub_repartition(CellMetric m, CellCellMetric ccm) override;

    //
    // Additional methods
    //

    // Clears obsolete entries from "partition"
    void clear_unknown_cell_ownership();

    /** Computes vector of vectors of cells to be sent to neighboring processes
     * based on the send volume passed in "send_volume".
     */
    PerNeighbor<GlobalCellIndices>
    compute_send_list(std::vector<double> &&send_volume,
                      const std::vector<double> &weights) const;

    /** Generate neighborhood information to send to processes that received
     * the cells in "cells_to_send".
     */
    PerNeighbor<__diff_impl::CellNeighborhoodPerCell>
    get_neighborhood_information(
        const PerNeighbor<GlobalCellIndices> &cells_to_send) const;

    /** Update "partition" vector based on cell neighborhood information
     * from neighboring processes.
     */
    void update_partitioning_from_received_neighbourhood(
        const PerNeighbor<__diff_impl::CellNeighborhoodPerCell> &neighbourhood);

    const std::unordered_map<std::string, diff_variants::FlowCalcKind>
        supported_default_diffusion_variants
        = {{"willebeek", diff_variants::FlowCalcKind::WILLEBEEK},
           {"schornbaum", diff_variants::FlowCalcKind::SCHORN}};
};
} // namespace grids
} // namespace repa
