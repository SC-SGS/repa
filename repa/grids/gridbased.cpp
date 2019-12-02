
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

#include "gridbased.hpp"

#include <algorithm>
#include <boost/mpi/collectives.hpp>
#include <regex>

#include "util/ensure.hpp"
#include "util/mpi_cart.hpp"
#include "util/mpi_graph.hpp"
#include "util/push_back_unique.hpp"
#include "util/vdist.hpp"

#ifndef NDEBUG
#define GRID_DEBUG
#endif

namespace repa {
namespace grids {

rank_type GridBasedGrid::gloidx_to_rank(global_cell_index_type idx)
{
    return position_to_rank(gbox.midpoint(idx));
}

std::array<Vec3d, 8> GridBasedGrid::bounding_box(rank_type r)
{
    const Vec3i coord = util::mpi_cart_get_coords(comm_cart, r);
    const Vec3i dims = util::mpi_cart_get_dims(comm_cart);

    std::array<Vec3d, 8> result;
    size_t i = 0;
    // Ranks holding the bounding box grid points of "r" = (c0, c1, c2) are:
    // (c0,     c1,     c2) upper right back corner,
    // (c0 - 1, c1,     c2) upper left back corner,
    // (c0,     c1 - 1, c2) lower right back corner,
    // (c0,     c1,     c2 - 1) upper right front corner,
    // (c0 - 1, c1 - 1, c2) lower left back corner
    // ... 2 more ...
    // (c0 - 1, c1 - 1, c2 - 1) lower left front corner
    // In total the set: {c0, c0 - 1} x {c1, c1 - 1} x {c2, c2 - 1}
    Vec3i off;
    for (off[0] = 0; off[0] <= 1; ++off[0]) {
        for (off[1] = 0; off[1] <= 1; ++off[1]) {
            for (off[2] = 0; off[2] <= 1; ++off[2]) {
                Vec3i nc, mirror{0, 0, 0};

                for (int d = 0; d < 3; ++d) {
                    nc[d] = coord[d] - off[d];

                    // Periodically wrap to the correct processor
                    // and save the wrapping to correct the grid point later.
                    // Can only happen in negative direction.
                    if (nc[d] < 0) {
                        nc[d] = dims[d] - 1;
                        mirror[d] = -1;
                    }
                }

                rank_type proc = util::mpi_cart_rank(comm_cart, nc);

                // Mirror the gridpoint back to where this subdomain is
                // expecting it.
                for (int d = 0; d < 3; ++d)
                    result[i][d] = gridpoints[proc][d] + mirror[d] * box_l[d];
                i++;
            }
        }
    }
    return result;
}

void GridBasedGrid::init_partitioning()
{
    is_regular_grid = true;

    // Copy data from grid.hpp
    for (int d = 0; d < 3; ++d)
        gridpoint[d] = (node_pos[d] + 1) * (box_l[d] / node_grid[d]);

    init_neighbors();
    init_octagons();
}

void GridBasedGrid::init_neighbors()
{
    neighbor_ranks.clear();
    neighbor_idx.clear();

    const Vec3i coord = util::mpi_cart_get_coords(comm_cart);
    const Vec3i dims = util::mpi_cart_get_dims(comm_cart);

    std::vector<rank_type> source_neigh,
        dest_neigh; // Send and receive neighborhood for repart
    rank_index_type nneigh = 0;
    Vec3i off;
    for (off[0] = -1; off[0] <= 1; ++off[0]) {
        for (off[1] = -1; off[1] <= 1; ++off[1]) {
            for (off[2] = -1; off[2] <= 1; ++off[2]) {
                Vec3i nc;

                for (int d = 0; d < 3; ++d) {
                    nc[d] = coord[d] + off[d];

                    // Periodic wrap
                    if (nc[d] < 0)
                        nc[d] = dims[d] - 1;
                    else if (nc[d] == dims[d])
                        nc[d] = 0;
                }

                rank_type r = util::mpi_cart_rank(comm_cart, nc);

                // Insert "r" as a new neighbor if yet unseen.
                if (r == comm_cart.rank())
                    continue;
                if (std::find(std::begin(neighbor_ranks),
                              std::end(neighbor_ranks), r)
                    == std::end(neighbor_ranks)) {
                    neighbor_ranks.push_back(r);
                    nneigh++;
                }

                if (off[0] >= 0 && off[1] >= 0 && off[2] >= 0)
                    util::push_back_unique(source_neigh, r);
                if (off[0] <= 0 && off[1] <= 0 && off[2] <= 0)
                    util::push_back_unique(dest_neigh, r);
            }
        }
    }

    std::sort(std::begin(neighbor_ranks), std::end(neighbor_ranks));
    // Inverse mapping
    for (rank_index_type i = 0; i < nneigh; ++i)
        neighbor_idx[neighbor_ranks[i]] = i;

    source_neigh.push_back(comm_cart.rank());
    dest_neigh.push_back(comm_cart.rank());
    neighcomm
        = util::directed_graph_communicator(comm_cart, source_neigh, dest_neigh);
}

void GridBasedGrid::init_octagons()
{
    boost::mpi::all_gather(comm_cart, gridpoint, gridpoints);

    my_dom = util::tetra::Octagon(bounding_box(comm_cart.rank()));

    neighbor_doms.clear();
    neighbor_doms.reserve(neighbor_ranks.size());

    for (size_t i = 0; i < neighbor_ranks.size(); ++i) {
        neighbor_doms.push_back(
            util::tetra::Octagon(bounding_box(neighbor_ranks[i])));
    }
}

void GridBasedGrid::reinit()
{
    nlocalcells = 0;
    nghostcells = 0;
    cells.clear();
    global_to_local.clear();
    exchange_vec.clear();

    // Reinit cells, nlocalcells, global_to_local
    // Simple loop over all global cells; TODO: optimize
    for (global_cell_index_type i = 0; i < gbox.ncells(); ++i) {
        auto midpoint = gbox.midpoint(i);

        // We use "position_to_rank" here.
        // Note that .contains() can be true for several octagons (see
        // comment in "position_to_rank"). Hence, we use position_to_rank
        // as tie-breaker.
        // We, however, use my_dom.contains() here as a guard so
        // "position_to_rank" does not throw.
        if (my_dom.contains(midpoint)
            && position_to_rank(midpoint) == comm.rank()) {
            cells.push_back(i);
            global_to_local[i] = nlocalcells;
            nlocalcells++;
        }
    }

#ifdef GRID_DEBUG
    printf("[%i] nlocalcells: %i\n", comm_cart.rank(), nlocalcells);
#endif
    ENSURE(nlocalcells > 0);

    // Temporary storage for exchange descriptors.
    // Will be filled only for neighbors
    // and moved from later.
    exchange_vec.clear();
    exchange_vec.resize(neighbor_ranks.size());

    // Determine ghost cells and communication volume
    for (local_cell_index_type i = 0; i < nlocalcells; i++) {
        for (global_cell_index_type neighidx :
             gbox.full_shell_neigh_without_center(cells[i])) {
            rank_type owner = gloidx_to_rank(neighidx);

            if (owner == comm_cart.rank())
                continue;

            // Add ghost cells only once to "cells" vector.
            if (global_to_local.find(neighidx) == std::end(global_to_local)) {
                // Add ghost cell to cells vector
                cells.push_back(neighidx);
                // Index mapping from global to ghost
                global_to_local[neighidx] = nlocalcells + nghostcells;
                // Number of ghost cells
                nghostcells++;
            }

            rank_index_type idx = neighbor_idx[owner];
            // Initialize exdesc and add "rank" as neighbor if unknown.
            if (exchange_vec[idx].dest == -1)
                exchange_vec[idx].dest = owner;

            util::push_back_unique(exchange_vec[idx].recv, neighidx);
            util::push_back_unique(exchange_vec[idx].send, cells[i]);
        }
    }

#ifdef GRID_DEBUG
    printf("[%i] nghostcells: %i\n", comm_cart.rank(), nghostcells);
    ENSURE(comm_cart.size() == 1 || nghostcells > 0);
#endif

    // All neighbors must be communicated with, otherwise something went wrong.
    // Sort and global_to_local.
    const auto glo_to_loc = [this](global_cell_index_type i) {
#ifdef GRID_DEBUG
        return global_to_local.at(i);
#else
        return global_to_local[i];
#endif
    };

    // Note:
    // Due to the definition of cell-ownership it can happen, that a
    // exchange vector is empty.
    // Therefore, we delete the empty ones here.
    exchange_vec.erase(
        std::remove_if(std::begin(exchange_vec), std::end(exchange_vec),
                       [](const GhostExchangeDesc &r) { return r.dest == -1; }),
        std::end(exchange_vec));

    for (auto &v : exchange_vec) {
        ENSURE(v.dest != -1);

        std::sort(std::begin(v.recv), std::end(v.recv));
        std::transform(std::begin(v.recv), std::end(v.recv), std::begin(v.recv),
                       glo_to_loc);

        std::sort(std::begin(v.send), std::end(v.send));
        std::transform(std::begin(v.send), std::end(v.send), std::begin(v.send),
                       glo_to_loc);
    }
}

GridBasedGrid::GridBasedGrid(const boost::mpi::communicator &comm,
                             Vec3d box_size,
                             double min_cell_size,
                             ExtraParams ep)
    : ParallelLCGrid(comm, box_size, min_cell_size),
      mu(1.0),
      gbox(box_size, min_cell_size),
      subdomain_midpoint(ep.subdomain_midpoint
                             ? ep.subdomain_midpoint
                             : decltype(subdomain_midpoint){std::bind(
                                   &GridBasedGrid::get_subdomain_center, this)})
{
    init_partitioning();
    reinit();
}

GridBasedGrid::~GridBasedGrid()
{
}

local_cell_index_type GridBasedGrid::n_local_cells()
{
    return nlocalcells;
}

ghost_cell_index_type GridBasedGrid::n_ghost_cells()
{
    return nghostcells;
}

rank_index_type GridBasedGrid::n_neighbors()
{
    return neighbor_ranks.size();
}

rank_type GridBasedGrid::neighbor_rank(rank_index_type i)
{
    return neighbor_ranks[i];
}

local_or_ghost_cell_index_type
GridBasedGrid::cell_neighbor_index(local_cell_index_type cellidx,
                                   fs_neighidx neigh)
{
    return global_to_local.at(gbox.neighbor(cells[cellidx], neigh));
}

std::vector<GhostExchangeDesc> GridBasedGrid::get_boundary_info()
{
    return exchange_vec;
}

local_cell_index_type GridBasedGrid::position_to_cell_index(Vec3d pos)
{
    local_cell_index_type c;
    try {
        c = global_to_local.at(gbox.cell_at_pos(pos));
    }
    catch (...) {
        throw std::domain_error("Position not in local subdomain.");
    }
    // global_to_local also resolves ghost indices.
    if (c >= n_local_cells())
        throw std::domain_error("Position not in local subdomain.");
    return c;
}

rank_type GridBasedGrid::cart_topology_position_to_rank(Vec3d pos)
{
    Vec3i grid_coord;
    for (size_t i = 0; i < 3; ++i) {
        grid_coord[i] = pos[i] / (box_l[i] / node_grid[i]);
        if (grid_coord[i] >= node_grid[i])
            grid_coord[i] = node_grid[i] - 1;
        else if (grid_coord[i] < 0)
            grid_coord[i] = 0;
    }

    return util::mpi_cart_rank(comm_cart, grid_coord);
}

rank_type GridBasedGrid::position_to_rank(Vec3d pos)
{
    // Cell ownership is based on the cell midpoint.
    auto mp = gbox.midpoint(gbox.cell_at_pos(pos));

    // This is fragile. Hopefully, MPI_Cart_rank defines cell ownership
    // as we do below...
    if (is_regular_grid)
        return cart_topology_position_to_rank(mp);

    // Cell ownerership is not uniquely determined by .contains()
    // because this function also accepts points on the boundary of an octagon.
    //
    // We simply define the owner to be the one with the lowest rank
    // among all processes where .contains() evaluates to true.
    //
    // Note, that neighbor_ranks is ordered by rank.
    rank_index_type i;
    for (i = 0; i < n_neighbors() && neighbor_ranks[i] < comm.rank(); ++i) {
        if (neighbor_doms[i].contains(mp))
            return neighbor_ranks[i];
    }

    if (my_dom.contains(mp))
        return comm.rank();

    for (; i < n_neighbors(); ++i) {
        if (neighbor_doms[i].contains(mp))
            return neighbor_ranks[i];
    }

    throw std::domain_error("Position unknown. Possibly a position outside of "
                            "the neighborhood of this process.");
}

rank_index_type GridBasedGrid::position_to_neighidx(Vec3d pos)
{
    rank_type rank = position_to_rank(pos);
    try {
        return neighbor_idx.at(rank);
    }
    catch (...) {
        throw std::domain_error("Position not in neighboring subdomain.");
    }
}

Vec3d GridBasedGrid::cell_size()
{
    return gbox.cell_size();
}

Vec3i GridBasedGrid::grid_size()
{
    return gbox.grid_size();
}

Vec3d GridBasedGrid::get_subdomain_center()
{
    Vec3d c{0., 0., 0.};

    // If no particles: Use subdomain midpoint.
    // (Calculated as mispoint of all cells).
    for (local_cell_index_type i = 0; i < n_local_cells(); ++i) {
        auto mp = gbox.midpoint(cells[i]);
        for (int d = 0; d < 3; ++d)
            c[d] += mp[d];
    }

    const local_cell_index_type n = n_local_cells();
    for (int d = 0; d < 3; ++d)
        c[d] /= n;

    return c;
}

bool GridBasedGrid::repartition(CellMetric m,
                                CellCellMetric ccm,
                                Thunk exchange_start_callback)
{
    // The node displacement is calculated according to
    // C. Begau, G. Sutmann, Comp. Phys. Comm. 190 (2015), p. 51 - 61
    rank_index_type nneigh = util::mpi_undirected_neighbor_count(neighcomm);

    auto weights = m();
    if (weights.size() != n_local_cells()) {
        throw std::runtime_error(
            "Metric only supplied " + std::to_string(weights.size())
            + "weights. Necessary: " + std::to_string(n_local_cells()));
    }
    double lambda_p
        = std::accumulate(std::begin(weights), std::end(weights), 0.0);
    auto r_p = this->subdomain_midpoint();

    std::vector<double> lambda(nneigh);
    MPI_Neighbor_allgather(&lambda_p, 1, MPI_DOUBLE, lambda.data(), 1,
                           MPI_DOUBLE, neighcomm);

    double lnormalizer
        = std::accumulate(lambda.begin(), lambda.end(), 0.0) / nneigh;

    std::vector<double> lambda_hat(nneigh);
    for (int i = 0; i < nneigh; ++i)
        lambda_hat[i] = lambda[i] / lnormalizer;

    std::vector<double> r(3 * nneigh);
    MPI_Neighbor_allgather(r_p.data(), 3, MPI_DOUBLE, r.data(), 3, MPI_DOUBLE,
                           neighcomm);

    for (rank_index_type i = 0; i < nneigh; ++i) {
        // Form "u"
        for (int d = 0; d < 3; ++d)
            r[3 * i + d] -= gridpoint[d];
        double len = util::norm2(&r[3 * i]);

        // Form "f"
        for (int d = 0; d < 3; ++d)
            r[3 * i + d] = (lambda_hat[i] - 1) * r[3 * i + d] / len;
    }

    const Vec3i coords = util::mpi_cart_get_coords(comm_cart);
    const Vec3i dims = util::mpi_cart_get_dims(comm_cart);

    Vec3d new_c = gridpoint;
    for (int d = 0; d < 3; ++d) {
        // Shift only non-boundary coordinates
        if (coords[d] == dims[d] - 1)
            continue;
        for (rank_index_type i = 0; i < nneigh; ++i)
            new_c[d] += mu * r[3 * i + d];
    }

    // Note: Since we do not shift gridpoints over periodic boundaries,
    // f values from periodic neighbors are not considered.
    // (See if condition in above loop.)
    // Therefore, they do not need periodic mirroring.

    // Note 2: We do not need to consider neighbors multiple times even
    // if two processes neighbor themselves along multiple boundaries.
    // We have a Cartesian grid. That means that if a process
    // appears twice in the neighborhood, all do.
    // So we can safely neglect multiple neighbors.

#ifdef GRID_DEBUG
    std::cout << "[" << comm_cart.rank() << "] Old c: " << gridpoint[0] << ","
              << gridpoint[1] << "," << gridpoint[2] << std::endl;
    std::cout << "[" << comm_cart.rank() << "] New c: " << new_c[0] << ","
              << new_c[1] << "," << new_c[2] << std::endl;
#endif

    // Update gridpoint and gridpoints
    // Currently allgather. Can be done in 64 process neighborhood.
    gridpoint = new_c;

    auto old_gridpoints = gridpoints;
    gridpoints.clear();
    boost::mpi::all_gather(comm_cart, gridpoint, gridpoints);
    ENSURE(gridpoints.size() == comm_cart.size());

    // Check for admissibility of new grid.
    // We do not constrain the grid cells to be convex.
    // But the bare minimum that we have to enforce is that grid points do
    // not collide with each other.

    const auto cs = cell_size();
    auto min_cell_size = std::min(std::min(cs[0], cs[1]), cs[2]);

    int nconflicts = 0;

    auto bb = bounding_box(comm_cart.rank());

    for (size_t i = 0; i < bb.size(); ++i)
        for (size_t j = i + 1; j < bb.size(); ++j)
            if (util::dist2(bb[i], bb[j]) < 2 * min_cell_size)
                nconflicts++;

    MPI_Allreduce(MPI_IN_PLACE, &nconflicts, 1, MPI_INT, MPI_SUM, comm_cart);

    if (nconflicts > 0) {
        std::cout << "Gridpoint update rejected because of node conflicts."
                  << std::endl;
        ENSURE(0);
        gridpoints = old_gridpoints;
        gridpoint = gridpoints[comm_cart.rank()];
        return false;
    }

    is_regular_grid = false;

    init_octagons();
    exchange_start_callback();
    reinit();

    return true;
}

void GridBasedGrid::command(std::string s)
{
    static const std::regex mure("\\s*mu\\s*=\\s*(\\d+\\.|\\.\\d+|\\d+.\\d+)");
    std::smatch m;

    if (std::regex_match(s, m, mure)) {
        mu = std::strtod(m[1].str().c_str(), NULL);
        if (comm_cart.rank() == 0)
            std::cout << "Setting mu = " << mu << std::endl;
    }
}

global_cell_index_type
GridBasedGrid::global_hash(local_or_ghost_cell_index_type cellidx)
{
    return cells[cellidx];
}

} // namespace grids
} // namespace repa