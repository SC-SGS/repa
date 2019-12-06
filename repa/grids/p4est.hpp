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

#include "util/p4est_deleter.hpp"

#include "../pargrid.hpp"
#include <array>
#include <memory>
#include <vector>

namespace repa {
namespace grids {

namespace impl {

enum class CellType { inner = 0, boundary = 1, ghost = 2 };
struct CellInfo {
    const rank_type
        owner_rank;     // the rank of this cell (equals this_node for locals)
    CellType cell_type; // shell information (0: inner local cell, 1: boundary
                        // local cell, 2: ghost cell)
                        // boundary Cells with shell-type 2 are set if the are
                        // in the periodic halo
    std::array<int, 26>
        neighbor; // unique index of the fullshell neighborhood cells (as in
                  // p8est); only 26 because cell itself is not included.
    const Vec3i coord; // cartesian coordinates of the cell
                       // For a globally unique index (on a virtual regular
                       // grid) call impl::cell_morton_index(coord).

    CellInfo(rank_type rank, CellType shell, const Vec3i &coord)
        : owner_rank(rank), cell_type(shell), coord(coord)
    {
        neighbor.fill(-1);
    }
};

struct RepartState {
    const boost::mpi::communicator &comm_cart;
    bool after_repart;
    std::vector<p4est_locidx_t> nquads_per_proc;
    std::function<void()> exchange_start_callback;

    RepartState(const boost::mpi::communicator &comm_cart);

    inline void reset();
    inline void inc_nquads(rank_type proc);
    inline void allreduce();
};

} // namespace impl

struct P4estGrid : public ParallelLCGrid {
    P4estGrid(const boost::mpi::communicator &comm,
              Vec3d box_size,
              double min_cell_size);
    local_cell_index_type n_local_cells() override;
    ghost_cell_index_type n_ghost_cells() override;
    rank_index_type n_neighbors() override;
    rank_type neighbor_rank(rank_index_type i) override;
    Vec3d cell_size() override;
    Vec3i grid_size() override;
    local_or_ghost_cell_index_type
    cell_neighbor_index(local_cell_index_type cellidx,
                        fs_neighidx neigh) override;
    std::vector<GhostExchangeDesc> get_boundary_info() override;
    local_cell_index_type position_to_cell_index(Vec3d pos) override;
    rank_type position_to_rank(Vec3d pos) override;
    rank_index_type position_to_neighidx(Vec3d pos) override;
    bool repartition(CellMetric m,
                     CellCellMetric ccm,
                     Thunk exchange_start_callback) override;
    global_cell_index_type
    global_hash(local_or_ghost_cell_index_type cellidx) override;

private:
    int m_grid_level;

    // Number of grid cells in total and per tree.
    Vec3i m_grid_size, m_brick_size;
    // Cell size (box_l / m_grid_size)
    Vec3d m_cell_size, m_inv_cell_size;

    // p4est data structures
    std::unique_ptr<p8est_connectivity_t> m_p8est_conn;
    std::unique_ptr<p8est_t> m_p8est;
    local_cell_index_type m_num_local_cells;
    ghost_cell_index_type m_num_ghost_cells;

    // helper data structures
    std::vector<global_cell_index_type>
        m_node_first_cell_idx; // First morton index on a virtual regular grid
                               // spanning all processes.
#ifdef GLOBAL_HASH_NEEDED
    std::vector<global_cell_index_type>
        m_global_idx; //< Global virtual morton index of cells (for
                      // global_hash()) could be removed in non-test builds.
#endif
    std::vector<impl::CellInfo> m_p8est_cell_info; //< Raw cell info from p4est

    // comm data structures
    std::vector<GhostExchangeDesc> m_exdescs;
    std::vector<rank_type> m_neighranks;

    void set_optimal_cellsize();
    void create_grid();
    void prepare_communication();
    // Reinitialized the grid (instantiation or after repartitioning)
    void reinitialize();

    impl::RepartState m_repartstate;
};
} // namespace grids
} // namespace repa
