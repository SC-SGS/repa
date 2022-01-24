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
#include "util/global_index_storage.hpp"
#include "util/initial_partitioning.hpp"
#include "util/ioptional.hpp"
#include <unordered_map>

namespace repa {
namespace grids {

struct GloMethod : public ParallelLCGrid {
    GloMethod(const boost::mpi::communicator &comm,
              Vec3d box_size,
              double min_cell_size,
              ExtraParams ep);
    virtual ~GloMethod();
    void after_construction() override;
    util::const_span<rank_type> neighbor_ranks() const override;
    Vec3d cell_size() const override;
    Vec3i grid_size() const override;
    local_or_ghost_cell_index_type
    cell_neighbor_index(local_cell_index_type cellidx,
                        fs_neighidx neigh) override;
    util::const_span<GhostExchangeDesc> get_boundary_info() override;
    local_cell_index_type position_to_cell_index(Vec3d pos) override;
    rank_type position_to_rank(Vec3d pos) override;
    bool repartition(CellMetric m,
                     CellCellMetric ccm,
                     Thunk exchange_start_callback) override;

    global_cell_index_type
    global_hash(local_or_ghost_cell_index_type cellidx) override;

protected:
    friend class HybridGPDiff;

    local_cell_index_type n_local_cells() const override;
    ghost_cell_index_type n_ghost_cells() const override;

    virtual bool sub_repartition(CellMetric m, CellCellMetric ccm) = 0;

    /** All neighbor ranks (ranks of subdomains neighboring this subdomain)
     * Do not use raw access to this vector in the implementation of GloMethod.
     * Use n_neighbors() and neighbor_rank(rank_index_type) instead.
     * At least one subclass (GridBasedGrid) uses its own implementaion of the
     * aforementioned menthods which is inconsistent to the ones stored in this
     * vector (might be sorted differently *and* "neighbors" is not available
     * during initialization).
     */
    std::vector<rank_type> neighbors;

    // Communication descriptors
    std::vector<GhostExchangeDesc> exchangeVector;

    // Defines the linearization and global cell neighborhoods.
    globox::GlobalBox<global_cell_index_type, int> gbox;

    // Stores the initial partitioning enum
    const util::InitialPartitionType initial_partitioning;

    // Stores the global index of local cells and the ghost cells of the
    // subdomain associated to this process.
    util::global_index_storage cell_store;

    // Called before init()
    // To override by subclasses if necessary
    // Partitioning is done, however all datastructures (but "partition")
    // are still in the old state (inconsistent to the new partitioning).
    virtual void pre_init(bool firstcall)
    {
    }

    // Called after init()
    // To override by subclasses if necessary
    // Partitioning and init fully done.
    virtual void post_init(bool firstcall)
    {
    }

    // Callback for a new border cell
    // Called during init()
    // @param localcell local cell index of the border cell
    // @apram foreigncell global index of the new ghost cell
    // @param owner owner rank of the ghost cell
    virtual void init_new_foreign_cell(local_cell_index_type localcell,
                                       global_cell_index_type foreigncell,
                                       rank_type owner)
    {
    }

    // For re-initialization, computes the new set of local cells
    // owned by this process/subdomain.
    // Override this if the subclass can do this more efficiently.
    // @returns Vector of global indices of the new local cells
    virtual std::vector<global_cell_index_type> compute_new_local_cells() const;

    // Function that globally resolves a cell index to a rank.
    // Return value may be empty if cell "idx" is irrelevant for
    // this process (no local cell and no ghost cell).
    //
    //
    // @param idx cell index to resolve
    // @returns rank which is responsible for global cell index idx
    virtual util::ioptional<rank_type>
    rank_of_cell(global_cell_index_type idx) const = 0;

    // Reinitializes the subdomain and communication data structures
    // after repartitioning.
    void init(bool firstcall = false);
};
} // namespace grids
} // namespace repa
