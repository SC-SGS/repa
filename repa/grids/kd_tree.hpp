/**
 * Copyright 2017-2019 Steffen Hirschmann
 * Copyright 2017-2018 Adriaan Nieß
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
#include <utility>
#include <vector>

#include "_compat.hpp"
#include "pargrid.hpp"

namespace repa {
namespace grids {

/**
 * A domain is a 3d-box defined by a tuple of vectors for the lower corner
 * and upper corner. While the domain includes the coordinate of the lower
 * corner, it excludes the coordinate of the upper corner.
 */
using Domain = std::pair<Vec3i, Vec3i>;

/** Encapsulates a kdpart::PartTreeStorage from libkdpart in kd_tree.cpp
 * to not having to include kdpart.h in this header.
 */
struct KDTreePrivateImpl;

class KDTreeGrid : public ParallelLCGrid {
private:
    /** Size of the global simulation box in cells. */
    const Vec3i m_global_domain_size;

    /** Domain representing the global simulation box. */
    const Domain m_global_domain;

    const Domain m_global_ghostdomain;
    const Vec3d m_cell_size;

    /** Internal k-d tree datastructure. */
    std::unique_ptr<KDTreePrivateImpl> m_kdtree;

    Domain m_local_subdomain;
    Domain m_local_ghostdomain;
    Vec3i m_local_subdomain_size;
    Vec3i m_local_ghostdomain_size;
    local_cell_index_type m_nb_of_local_cells;
    ghost_cell_index_type m_nb_of_ghost_cells;
    std::vector<local_or_ghost_cell_index_type> m_index_permutations;
    std::vector<local_or_ghost_cell_index_type> m_index_permutations_inverse;

    /** Maps neighbor id (nidx) to rank. */
    std::vector<rank_type> m_neighbor_processes;

    /** Maps rank to neighbor id (nidx) or -1 if rank is no neighbor. */
    std::vector<rank_index_type> m_neighbor_processes_inverse;

    std::vector<GhostExchangeDesc> m_boundary_info;

    void init_local_domain_bounds();
    void init_nb_of_cells();
    void init_index_permutations();

    /** Transforms domain vector to cellvector. */
    std::vector<Vec3i> cells(const std::vector<Domain> &domains);

    /**
     * Initializes datastructure that contains the ranks of all neighbor
     * processes.
     * TODO update comment
     */
    void init_neighborhood_information();

    void init_neighborhood_information(int neighbor_rank);

    /** Init receiving ghostcells. */
    void init_recv_cells(GhostExchangeDesc &gexd,
                         const Domain &neighbor_subdomain);

    /** Init sending local cells. */
    void init_send_cells(GhostExchangeDesc &gexd,
                         const Domain &neighbor_ghostdomain);

    void clear_lookup_datastructures();

    void reinitialize();

public:
    KDTreeGrid(const boost::mpi::communicator &comm,
               Vec3d box_size,
               double min_cell_size);

    virtual local_cell_index_type n_local_cells() override;

    virtual ghost_cell_index_type n_ghost_cells() override;

    virtual rank_index_type n_neighbors() override;

    virtual rank_type neighbor_rank(rank_index_type i) override;

    virtual Vec3d cell_size() override;

    virtual Vec3i grid_size() override;

    virtual local_or_ghost_cell_index_type
    cell_neighbor_index(local_cell_index_type cellidx,
                        fs_neighidx neigh) override;

    virtual std::vector<GhostExchangeDesc> get_boundary_info() override;

    virtual local_cell_index_type position_to_cell_index(Vec3d pos) override;

    virtual rank_type position_to_rank(Vec3d pos) override;

    /**
     * @throws std::domain_error if position is not on a neighboring process.
     */
    virtual rank_index_type position_to_neighidx(Vec3d pos) override;

    virtual bool
    repartition(CellMetric m, CellCellMetric ccm, Thunk cb) override;

    global_cell_index_type
    global_hash(local_or_ghost_cell_index_type cellidx) override;
};

} // namespace grids
} // namespace repa
