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

#include "glomethod.hpp"
#include "util/tetra.hpp"

namespace repa {
namespace grids {

/** Implements a grid-based load-balancing scheme.
 * A regular grid partitioning grid is layed over the regular cell grid.
 * The vertices of the partitioning grid are shifted towards local load
 * centers for overloaded subdomains.
 * This keeps the communication structure between processes constant.
 */
struct GridBasedGrid : public GloMethod {
    GridBasedGrid(const boost::mpi::communicator &comm,
                  Vec3d box_size,
                  double min_cell_size,
                  ExtraParams ep);
    ~GridBasedGrid() override;
    void command(std::string s) override;

    // Neighborhood is not determined via GloMethod. Because of the constancy
    // (see below), this grid implementation provides its own implementation of
    // neighborhood. DO NOT REMOVE these. GridBasedGrid::rank_of_cell needs this
    // neighborhood.
    rank_index_type n_neighbors() const override;
    rank_type neighbor_rank(rank_index_type i) const override;

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

    // Triangulation data structure for this subdomain
    util::tetra::Octagon my_dom;

    /** Note that the number of neighbors
     * is constant and the neighbors themselves do *not*
     * change over time.
     * However, since we do not know how many neighbors a subdomain
     * will have (i.e. if nproc < 26 vs. nproc is prime vs. nproc = 10^3)
     * we still do not keep them in a rigid data structure of size 26.
     */
    std::vector<rank_type> const_neighborhood;
    // Triangulation data structure for the neighboring subdomains
    std::vector<util::tetra::Octagon> neighbor_doms;

    // Associated grid point -- upper right back vertex of subdomain.
    Vec3d gridpoint;
    // The gathered version of "gridpoint", i.e. the gridpoint of every process.
    std::vector<Vec3d> gridpoints;

    // Neighborhood communicator for load exchange during repart
    boost::mpi::communicator neighcomm;

    // Returns the 8 vertices bounding the subdomain of rank "r"
    std::array<Vec3d, 8> bounding_box(rank_type r) const;

    bool sub_repartition(CellMetric m, CellCellMetric ccm) override;
    util::ioptional<rank_type>
    rank_of_cell(global_cell_index_type idx) const override;
    void pre_init(bool firstcall) override;
    void post_init(bool firstcall) override;

    // Initializes the partitioning to a regular Cartesian grid.
    void init_regular_partitioning();

    // Initializes "my_dom" and "neighbor_doms"
    void init_octagons();

    // Initializes the neighbor ranks data structures
    void create_cartesian_neighborhood();

    // Returns the center of this subdomain
    Vec3d get_subdomain_center();

    // Check if shifted gridpoint is valid
    bool check_validity_of_subdomains(const std::vector<rank_type> &) const;

    // Resolves a position on the initial Cartesian decomposition, set
    // by init_regular_partitioning.
    rank_type cart_topology_position_to_rank(Vec3d pos) const;

    // Function returning the subdomain midpoint.
    // Either a user-passed function via ExtraParams in the constructor
    // or get_subdomain_center().
    decltype(ExtraParams::subdomain_midpoint) subdomain_midpoint;
};
} // namespace grids
} // namespace repa
