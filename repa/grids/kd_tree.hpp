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
#include "globox.hpp"
#include "grid_variants.hpp"
#include "pargrid.hpp"
#include "util/box_global_index_storage.hpp"

namespace repa {
namespace grids {

/** Encapsulates a kdpart::PartTreeStorage from libkdpart in kd_tree.cpp
 * to not having to include kdpart.h in this header.
 */
struct KDTreePrivateImpl;

class KDTreeGrid : public ParallelLCGrid, public VariantSetter {
public:
    KDTreeGrid(const boost::mpi::communicator &comm,
               Vec3d box_size,
               double min_cell_size);
    virtual util::const_span<rank_type> neighbor_ranks() const override;
    virtual Vec3d cell_size() const override;
    virtual Vec3i grid_size() const override;
    virtual local_or_ghost_cell_index_type
    cell_neighbor_index(local_cell_index_type cellidx,
                        fs_neighidx neigh) override;
    virtual util::const_span<GhostExchangeDesc> get_boundary_info() override;
    virtual local_cell_index_type position_to_cell_index(Vec3d pos) override;
    virtual rank_type position_to_rank(Vec3d pos) override;
    virtual bool
    repartition(CellMetric m, CellCellMetric ccm, Thunk cb) override;
    global_cell_index_type
    global_hash(local_or_ghost_cell_index_type cellidx) override;

    virtual void command(std::string s);

    virtual std::set<std::string> get_supported_variants() const override;
    virtual void set_variant(const std::string &var) override;

    /**
     * A domain is a 3d-box defined by a tuple of vectors for the lower corner
     * and upper corner. While the domain includes the coordinate of the lower
     * corner, it excludes the coordinate of the upper corner.
     */
    using Domain = std::pair<Vec3i, Vec3i>;

protected:
    virtual local_cell_index_type n_local_cells() const override;
    virtual ghost_cell_index_type n_ghost_cells() const override;

private:
    /** Internal k-d tree datastructure. */
    std::unique_ptr<KDTreePrivateImpl> m_kdtree;

    Domain m_local_subdomain;
    Domain m_local_ghostdomain;
    Vec3i m_local_ghostdomain_size;

    util::box_global_index_storage cell_store;

    /** Maps neighbor id (nidx) to rank. */
    std::vector<rank_type> m_neighbor_processes;

    globox::GlobalBox<global_cell_index_type, int> gbox;
    std::vector<GhostExchangeDesc> m_boundary_info;

    /** Flag whether to perform local or global repartitioning,
     * modified via command().
     */
    bool local_repart = false;

    /**
     * Initializes datastructure that contains neighbor ranks and ghost
     * communication
     */
    void init_neighborhood_information();

    void init_neighborhood_information(int neighbor_rank);

    /** Init receiving ghostcells. */
    void init_recv_cells(GhostExchangeDesc &gexd,
                         const Domain &neighbor_subdomain);

    /** Init sending local cells. */
    void init_send_cells(GhostExchangeDesc &gexd,
                         const Domain &neighbor_ghostdomain);

    void reinitialize();
};

} // namespace grids
} // namespace repa
