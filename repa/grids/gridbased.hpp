/**
 * Copyright 2017-2019 Steffen Hirschmann
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

#include <mpi.h>
#include <unordered_map>
#include <vector>

#include "globox.hpp"
#include "pargrid.hpp"
#include "util/tetra.hpp"

namespace repa {
namespace grids {

/** Implements a grid-based load-balancing scheme.
 * A regular grid partitioning grid is layed over the regular cell grid.
 * The vertices of the partitioning grid are shifted towards local load
 * centers for overloaded subdomains.
 * This keeps the communication structure between processes constant.
 */
struct GridBasedGrid : public ParallelLCGrid {
    GridBasedGrid(const boost::mpi::communicator &comm,
                  Vec3d box_size,
                  double min_cell_size,
                  ExtraParams ep);
    ~GridBasedGrid();
    void after_construction() override;
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
    /**
     * @throws std::domain_error if "pos" is not on a neighboring process.
     */
    rank_type position_to_rank(Vec3d pos) override;
    rank_index_type position_to_neighidx(Vec3d pos) override;
    bool repartition(CellMetric m,
                     CellCellMetric ccm,
                     Thunk exchange_start_callback) override;

    void command(std::string s) override;

    global_cell_index_type
    global_hash(local_or_ghost_cell_index_type cellidx) override;

private:
    // Indicator if the decomposition currently is a regular grid,
    // which is the case directly after instantiation.
    // This is important for position-to-rank queries.
    // They can be answered for the whole domain if the grid
    // is a regular grid. Otherwise, a process can only resolve
    // positions in neighboring subdomains.
    bool is_regular_grid;

    // Factor for grid point displacement.
    // Settable via command()
    double mu;

    // Number of local and ghost cells
    local_cell_index_type nlocalcells;
    ghost_cell_index_type nghostcells;

    // Triangulation data structure for this subdomain
    util::tetra::Octagon my_dom;

    // Triangulation data structure for the neighboring subdomains
    std::vector<util::tetra::Octagon> neighbor_doms;
    // Ranks of the neigbors. Note that the number of neighbors
    // is constant and the neighbors themselves do *not*
    // change over time.
    // However, since we do not know how many neighbors a subdomain
    // will have (i.e. if nproc < 26 vs. nproc is prime vs. nproc = 10^3)
    // we use a dynamic std::vector here.
    std::vector<rank_type> neighbor_ranks;

    // Inverse mapping neighbor_rank to index [0, 26) in "neighbor_ranks".
    std::unordered_map<rank_type, rank_index_type> neighbor_idx;

    // Associated grid point -- upper right back vertex of subdomain.
    Vec3d gridpoint;
    // The gathered version of "gridpoint", i.e. the gridpoint of every process.
    std::vector<Vec3d> gridpoints;

    // Indices of locally known cells. Local cells before ghost cells.
    std::vector<global_cell_index_type> cells;

    globox::GlobalBox<global_cell_index_type, global_cell_index_type> gbox;
    // Global to local index mapping, defined for local and ghost cells
    std::unordered_map<global_cell_index_type, local_cell_index_type>
        global_to_local;

    std::vector<GhostExchangeDesc> exchange_vec;

    // Returns the 8 vertices bounding the subdomain of rank "r"
    std::array<Vec3d, 8> bounding_box(rank_type r);

    // Neighborhood communicator for load exchange during repart
    boost::mpi::communicator neighcomm;

    // Global cell index to rank mapping
    rank_type gloidx_to_rank(global_cell_index_type idx);

    // Initializes the partitioning to a regular Cartesian grid.
    void init_partitioning();
    // Checks is "pos" is also accepted by a neighboring octagon.
    bool does_neighbor_accept(Vec3d pos);
    // Reinitializes the internal data of this class
    void reinit();
    // Initializes "my_dom" and "neighbor_doms"
    void init_octagons();

    // Initializes the neighbor ranks data structures
    void init_neighbors();

    // Returns the center of this subdomain
    Vec3d get_subdomain_center();

    // Check if shifted gridpoint is valid
    bool check_validity_of_subdomains(const std::vector<rank_type> &);

    rank_type cart_topology_position_to_rank(Vec3d pos);

    // Function returning the contribution of a single cell to the subdomain
    // midpoint. Either a user-passed function via ExtraParams in the
    // constructor or as a default returns the midpoint of a cell
    decltype(ExtraParams::subdomain_center_contribution_of_cell)
        get_subdomain_center_contribution_of_cell;
};
} // namespace grids
} // namespace repa
