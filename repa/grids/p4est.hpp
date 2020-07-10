/**
 * Copyright 2017-2019 Steffen Hirschmann
 * Copyright 2017-2018 Malte Brunn
 * Copyright 2017-2018 Michael Lahnert
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

#include "../pargrid.hpp"
#include <memory>

namespace repa {
namespace grids {

struct P4estGrid : public ParallelLCGrid {
    P4estGrid(const boost::mpi::communicator &comm,
              Vec3d box_size,
              double min_cell_size);
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
    local_cell_index_type n_local_cells() const override;
    ghost_cell_index_type n_ghost_cells() const override;

private:
    struct _P4estGrid_impl;
    std::unique_ptr<_P4estGrid_impl> _impl;
};
} // namespace grids
} // namespace repa
