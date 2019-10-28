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

//#ifdef HAVE_P4EST

#include "util/p4est_deleter.hpp"

#include "../pargrid.hpp"
#include <array>
#include <memory>
#include <vector>

namespace repa {
namespace grids {

namespace impl {

enum class CellType { inner = 0, boundary = 1, ghost = 2 };
struct LocalShell {
    global_cell_index_type
        idx; // a unique index within all cells (as used by p8est for locals)
    rank_type which_proc; // the rank of this cell (equals this_node for locals)
    CellType shell; // shell information (0: inner local cell, 1: boundary local
                    // cell, 2: ghost cell)
    int boundary;   // Bit mask storing boundary info. MSB ...
                    // z_r,z_l,y_r,y_l,x_r,x_l LSB Cells with shell-type 0 or
                    // those located within the domain are always 0 Cells with
                    // shell-type 1 store information about which face is a
    // boundary Cells with shell-type 2 are set if the are in the
    // periodic halo
    std::array<int, 26>
        neighbor; // unique index of the fullshell neighborhood cells (as in
                  // p8est); only 26 because cell itself is not included.
    Vec3i coord;  // cartesian coordinates of the cell

    LocalShell(global_cell_index_type idx,
               rank_type rank,
               CellType shell,
               int boundary,
               int x,
               int y,
               int z)
        : idx(idx),
          which_proc(rank),
          shell(shell),
          boundary(boundary),
          coord(x, y, z)
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
    std::vector<global_cell_index_type> m_node_first_cell_idx;
    std::vector<impl::LocalShell> m_p8est_shell;
#ifdef GLOBAL_HASH_NEEDED
    std::vector<global_cell_index_type>
        m_global_idx; //< Global virtual morton index of cells (for
                      // global_hash()) could be removed in non-test builds.
#endif

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

//#endif // HAVE_P4EST