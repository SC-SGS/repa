
/**
 * Copyright 2017-2020 Steffen Hirschmann
 * Copyright 2020 Benjamin Vier
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

#include "util/mpi_cart.hpp"
#include "util/mpi_graph.hpp"
#include "util/push_back_unique.hpp"
#include "util/vdist.hpp"
#include "util/vec_arith.hpp"

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
    for (off[2] = 0; off[2] <= 1; ++off[2]) {
        for (off[1] = 0; off[1] <= 1; ++off[1]) {
            for (off[0] = 0; off[0] <= 1; ++off[0]) {
                using namespace util::vector_arithmetic;
                Vec3i nc = (coord - off) % dims;
                rank_type proc = util::mpi_cart_rank(comm_cart, nc);

                // Mirror the gridpoint back to where this subdomain is
                // expecting it.
                const Vec3i mirror
                    = -static_cast_vec<Vec3i>((coord == 0) && (off == 1));
                result[i] = gridpoints[proc] + mirror * box_l;
                i++;
            }
        }
    }
    return result;
}

void GridBasedGrid::init_partitioning()
{
    is_regular_grid = true;

    using namespace util::vector_arithmetic;
    gridpoint = (node_pos + 1) * (box_l / node_grid);

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
                using namespace util::vector_arithmetic;
                const Vec3i nc = (coord + off) % dims;
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

                if (all(off >= 0))
                    util::push_back_unique(source_neigh, r);
                if (all(off <= 0))
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
    neighcomm = util::directed_graph_communicator(comm_cart, source_neigh,
                                                  dest_neigh);
}

void GridBasedGrid::init_octagons()
{
    boost::mpi::all_gather(comm_cart, gridpoint, gridpoints);

    my_dom = util::tetra::Octagon(bounding_box(comm_cart.rank()));

    neighbor_doms.clear();
    neighbor_doms.reserve(neighbor_ranks.size());
    std::transform(std::begin(neighbor_ranks), std::end(neighbor_ranks),
                   std::back_inserter(neighbor_doms), [this](rank_type r) {
                       return util::tetra::Octagon(bounding_box(r));
                   });
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

        if (my_dom.contains(midpoint)) {
            cells.push_back(i);
            global_to_local[i] = nlocalcells;
            nlocalcells++;
        }
    }

#ifdef GRID_DEBUG
    printf("[%i] nlocalcells: %i\n", comm_cart.rank(), nlocalcells);
#endif
    assert(nlocalcells > 0);

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
    assert(comm_cart.size() == 1 || nghostcells > 0);
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
        assert(v.dest != -1);

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
    // Cache the octagons for each process.
    // They become invalid after the first repartitioning.
    // DO NOT CALL this function then.
    assert(is_regular_grid);

    static std::map<rank_type, util::tetra::Octagon> all_octs;
    for (rank_type i = 0; i < comm_cart.size(); ++i) {
        if (all_octs.find(i) == std::end(all_octs))
            all_octs[i] = util::tetra::Octagon(bounding_box(i));

        if (all_octs[i].contains(pos))
            return i;
    }

    throw std::domain_error(
        "Position globally unknown. This is a bug, please report it.");
}

rank_type GridBasedGrid::position_to_rank(Vec3d pos)
{
    // Cell ownership is based on the cell midpoint.
    const auto mp = gbox.midpoint(gbox.cell_at_pos(pos));

    // Directly invoked lambda expression to prohibit use of "pos" in future
    // changes. Resolving must be done for "mp".
    return [this](Vec3d mp) {
        // .contains() is mutually exclusive. The expectation is that most
        // queried positions belong to this node, so check it first. The order
        // of the neighbors is not relevant.
        if (my_dom.contains(mp))
            return comm.rank();

        for (rank_index_type i = 0; i < n_neighbors(); ++i) {
            if (neighbor_doms[i].contains(mp))
                return neighbor_ranks[i];
        }

        if (is_regular_grid)
            return cart_topology_position_to_rank(mp);

        throw std::domain_error(
            "Position unknown. Possibly a position outside of "
            "the neighborhood of this process.");
    }(mp);
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
    using namespace util::vector_arithmetic;
    Vec3d c{0., 0., 0.};

    for (local_cell_index_type i = 0; i < n_local_cells(); ++i)
        c += gbox.midpoint(cells[i]);
    return c / static_cast<double>(n_local_cells());
}

bool GridBasedGrid::repartition(CellMetric m,
                                CellCellMetric ccm,
                                Thunk exchange_start_callback)
{
    using namespace util::vector_arithmetic;

    // The node displacement is calculated according to
    // C. Begau, G. Sutmann, Comp. Phys. Comm. 190 (2015), p. 51 - 61
    const rank_index_type nneigh
        = util::mpi_undirected_neighbor_count(neighcomm);
    const auto weights = m();
    assert(weights.size() == n_local_cells());

    const double lambda_p
        = std::accumulate(std::begin(weights), std::end(weights), 0.0);
    const auto r_p = this->subdomain_midpoint();

    std::vector<double> lambda(nneigh);
    MPI_Neighbor_allgather(&lambda_p, 1, MPI_DOUBLE, lambda.data(), 1,
                           MPI_DOUBLE, neighcomm);

    const double lnormalizer
        = std::accumulate(lambda.begin(), lambda.end(), 0.0) / nneigh;

    std::vector<double> lambda_hat(nneigh);
    for (int i = 0; i < nneigh; ++i)
        lambda_hat[i] = lambda[i] / lnormalizer;

    std::vector<Vec3d> r(nneigh);
    MPI_Neighbor_allgather(r_p.data(), sizeof(Vec3d), MPI_BYTE, r.data(),
                           sizeof(Vec3d), MPI_BYTE, neighcomm);

    for (rank_index_type i = 0; i < nneigh; ++i) {
        // Form "u"
        r[i] -= gridpoint;
        const double len = util::norm2(r[i]);
        // Form "f"
        r[i] *= (lambda_hat[i] - 1) / len;
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

    // Update gridpoint and gridpoints
    // Currently allgather. Can be done in 64 process neighborhood.

    double factor = 1.;
    for (rank_type i = 0; i < 8; i++) {
        shift_gridpoint(r, factor, i);
        auto old_gridpoints = gridpoints;

        gridpoints.clear();
        boost::mpi::all_gather(comm_cart, gridpoint, gridpoints);
        assert(gridpoints.size() == comm_cart.size());

        // Check for admissibility of new grid.
        // Note: Don't change "my_dom" here in case we decide to reset to the
        // current state and return false (or, also reset if afterwards).
        const auto cs = cell_size();
        const auto max_cs = std::max(std::max(cs[0], cs[1]), cs[2]);
        int hasConflict
            = util::tetra::Octagon(bounding_box(comm_cart.rank()), max_cs)
                      .is_valid()
                  ? 0
                  : 1;
        MPI_Allreduce(MPI_IN_PLACE, &hasConflict, 1, MPI_INT, MPI_SUM,
                      comm_cart);

        if (hasConflict > 0) {
            std::cout << "Gridpoint update rejected because of node conflicts."
                      << std::endl;
            gridpoints = old_gridpoints;
            gridpoint = gridpoints[comm_cart.rank()];
            if (factor > 0.24) {
                i--;
                factor /= 2.;
            }
            else {
                return false;
            }
        }
        else {
            factor = 1.;
        }
    }
    is_regular_grid = false;

    init_octagons();
    exchange_start_callback();
    reinit();

    return true;
}

void GridBasedGrid::shift_gridpoint(std::vector<Vec3d> r,
                                   double factor,
                                   int iteration)
{
    const Vec3i coords = util::mpi_cart_get_coords(comm_cart);
    const Vec3i dims = util::mpi_cart_get_dims(comm_cart);
    const rank_type proc = util::mpi_cart_rank(comm_cart, coords);

    if (proc % 8 != iteration) {
        return;
    }

    auto old_gp = gridpoint;
    for (int d = 0; d < 3; ++d) {
        // Shift only non-boundary coordinates
        if (coords[d] == dims[d] - 1)
            continue;
        for (int i = 0; i < r.size(); i++)
            gridpoint[d] += factor * mu * r[i][d];
    }
#ifdef GRID_DEBUG
    std::cout << "[" << comm_cart.rank() << "] Old c: " << old_gp[0] << ","
              << old_gp[1] << "," << old_gp[2] << std::endl;
    std::cout << "[" << comm_cart.rank() << "] New c: " << gridpoint[0] << ","
              << gridpoint[1] << "," << gridpoint[2] << std::endl;
#endif
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