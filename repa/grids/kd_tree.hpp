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

//#ifdef HAVE_KDPART

#include <array>
#include <cstring>
#include <mpi.h>
#include <tuple>
#include <unordered_set>
#include <vector>

#include <kdpart.h>

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

class KDTreeGrid : public ParallelLCGrid {
private:
    /** Size of the global simulation box in cells. */
    const Vec3i m_global_domain_size;

    /** Domain representing the global simulation box. */
    const Domain m_global_domain;

    const Domain m_global_ghostdomain;
    const Vec3i m_global_ghostdomain_size;
    const Vec3d m_cell_dimensions;

    /** Internal k-d tree datastructure. */
    kdpart::PartTreeStorage m_kdtree;

    Domain m_local_subdomain;
    Domain m_local_ghostdomain;
    Vec3i m_local_subdomain_size;
    Vec3i m_local_ghostdomain_size;
    unsigned m_nb_of_local_cells;
    unsigned m_nb_of_ghost_cells;
    std::vector<lgidx> m_index_permutations;
    std::vector<lgidx> m_index_permutations_inverse;

    /** Maps neighbor id (nidx) to rank. */
    std::vector<int> m_neighbor_processes;

    /** Maps rank to neighbor id (nidx) or -1 if rank is no neighbor. */
    std::vector<int> m_neighbor_processes_inverse;

    std::vector<GhostExchangeDesc> m_boundary_info;

private:
    /** Returns the grid dimensions of the global simulation box in cells. */
    Vec3i grid_dimensions();

    /** Returns the cell size within the global simulation box. */
    Vec3d cell_dimensions(const Vec3i &grid_dimensions);

    /** Returns the number of cells from the size of a domain. */
    static int volume(Vec3i domain_size);

    /** Returns the number of cells from a given domain*/
    static int volume(Domain domain_bounds);

    /** Returns the ghostdomain from a given domain. */
    static Domain ghostdomain_bounds(const Domain &domain);

    /** Returns the size of a given domain */
    static Vec3i domain_size(const Domain &domain);

    /**
     * Returns true if the given cell is a ghostcell respective to the given
     * ghostdomain
     *
     * @param cell A cell vector relative to the given ghostdomain
     * @param ghostdomain A ghostdomain
     * @return True if cell is within the ghostlayer of the given ghostdomain
     */
    static bool is_ghost_cell(const Vec3i &cell, const Vec3i &ghostdomain_size);

    /** Returns true if the given domain contains the given cell vector. */
    static bool domain_contains_cell(const Domain &domain, const Vec3i &cell);

    /**
     * Transforms a global position within the the simulation box to a global
     * cell vector.
     */
    Vec3i absolute_position_to_cell_position(double absolute_position[3]);

    void init_local_domain_bounds();

    void init_nb_of_cells();

    void init_index_permutations();

    /**
     * This method returns the intersecting domains between a localdomain and
     * a ghostdomain. Multiple intersection domains are possible in case of a
     * periodic domain.
     *
     * @param localdomain A subdomain which doesn't exceed the bounds of the
     *  global domain.
     * @param ghostdomain A ghostdomain which doesn't exceed the global
     *  ghostdomain
     * @param ghostdomain_coords If this parameter value is true, then the
     *  lu-coordinate of the resulting intersection domains is relative to the
     *  ghostdomain parameter. Otherwise if this value is false, the lu-coord
     *  is relative to the bounds of the localdomain parameter.
     * @param periodic_intersections_only If this parameter is true, then only
     *  intersection domains are returned that are caused by the ghostdomain
     *  (provided by the ghostdomain parameter) exceeding the bounds of the
     *  domain along a periodic dimension. This parameter can be useful if e.g.
     *  the localdomain parameter is a subset of the ghostdomain parameter and
     *  only intersections caused by overlapping domains are of interest.
     * @return List of intersecting domains between the subdomain and the
     *  ghostdomain relative to the global domain.
     */
    std::vector<Domain> intersection_domains(const Domain &localdomain,
                                             const Domain &ghostdomain,
                                             bool ghostdomain_coords = false,
                                             bool periodic_intersections_only
                                             = false) const;

    /**
     * Returns true if the given localdomain and the given ghostdomain
     * intersect. This includes intersections that are the result of periodic
     * domain bounds.
     */
    bool are_domains_intersecting(const Domain &localdomain,
                                  const Domain &ghostdomain) const;

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

    virtual lidx n_local_cells() override;

    virtual gidx n_ghost_cells() override;

    virtual nidx n_neighbors() override;

    virtual rank neighbor_rank(nidx i) override;

    virtual Vec3d cell_size() override;

    virtual Vec3i grid_size() override;

    virtual lgidx cell_neighbor_index(lidx cellidx, int neigh) override;

    virtual std::vector<GhostExchangeDesc> get_boundary_info() override;

    virtual lidx position_to_cell_index(double pos[3]) override;

    virtual rank position_to_rank(double pos[3]) override;

    virtual nidx position_to_neighidx(double pos[3]) override;

    virtual bool
    repartition(CellMetric m, CellCellMetric ccm, Thunk cb) override;

    int global_hash(lgidx cellidx) override;
};

} // namespace grids
} // namespace repa

//#endif // HAVE_KDPART