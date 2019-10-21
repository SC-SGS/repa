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

#include <algorithm>
#include <boost/range/iterator_range.hpp>
#include <mpi.h>

#include "cart.hpp"
#include "util/linearize.hpp"
#include "util/neighbor_offsets.hpp"
#include "util/push_back_unique.hpp"
#include "util/vadd.hpp"

#include "_compat.hpp"

namespace repa {
namespace grids {

namespace impl {

static std::pair<Vec3i, Vec3i> determine_send_receive_bounds(
    const Vec3i &offset, int receive, const Vec3i &grid)
{
    Vec3i lc, hc;

    for (int i = 0; i < 3; ++i) {
        lc[i] = offset[i] <= 0 ? 1 : grid[i];
        hc[i] = offset[i] < 0 ? 1 : grid[i];

        // The receive area is actually in the ghost layer
        // so shift the corresponding indices.
        if (receive) {
            if (offset[i] > 0)
                lc[i] = hc[i] = grid[i] + 1;
            else if (offset[i] < 0)
                lc[i] = hc[i] = 0;
        }
    }
    return std::make_pair(lc, hc);
}

} // namespace impl

bool CartGrid::is_ghost_cell(const Vec3i &c)
{
    return c[0] == 0 || c[0] == m_ghost_grid_size[0] - 1 || c[1] == 0
           || c[1] == m_ghost_grid_size[1] - 1 || c[2] == 0
           || c[2] == m_ghost_grid_size[2] - 1;
}

rank CartGrid::proc_offset_to_rank(const Vec3i &offset)
{
    auto neighpos = util::vadd_mod(m_procgrid_pos, offset, m_procgrid);
    int rank;
    MPI_Cart_rank(comm_cart, neighpos.data(), &rank);
    return rank;
}

bool CartGrid::self_comm_necessary()
{
    return m_procgrid[0] < 2 || m_procgrid[1] < 2 || m_procgrid[2] < 2;
}

void CartGrid::fill_neighranks()
{
    m_neighranks.clear();

    for (const auto &offset : util::NeighborOffsets3D::raw) {
        // Push back unique neighbor ranks into m_neighbors
        auto rank = proc_offset_to_rank(offset);
        if (rank != comm.rank() || self_comm_necessary())
            util::push_back_unique(m_neighranks, rank);
    }
}

void CartGrid::create_index_permutations()
{
    int ncells = n_local_cells() + n_ghost_cells();
    m_to_pargrid_order.resize(ncells);
    m_from_pargrid_order.resize(ncells);

    int localidx = 0, ghostidx = n_local_cells();
    for (int i = 0; i < ncells; ++i) {
        auto c = util::unlinearize(i, m_ghost_grid_size);
        if (is_ghost_cell(c)) {
            m_from_pargrid_order[ghostidx] = i;
            m_to_pargrid_order[i] = ghostidx++;
        }
        else {
            m_from_pargrid_order[localidx] = i;
            m_to_pargrid_order[i] = localidx++;
        }
    }
}

void CartGrid::create_grid()
{
    for (int i = 0; i < 3; ++i) {
        // Copy infos from deprecated pargrid variables
        m_procgrid[i] = node_grid[i];
        m_procgrid_pos[i] = node_pos[i];

        // Local box info
        m_localbox[i] = box_l[i] / m_procgrid[i];
        m_lowerleft[i] = m_localbox[i] * m_procgrid_pos[i];

        // Grid and cell size
        if (max_range > ROUND_ERROR_PREC * box_l[0])
            m_grid_size[i] = static_cast<int>(m_localbox[i] / max_range);
        else
            m_grid_size[i] = 1;
        m_ghost_grid_size[i] = m_grid_size[i] + 2;

        m_cell_size[i] = m_localbox[i] / m_grid_size[i];
        m_inv_cell_size[i] = 1.0 / m_cell_size[i];
    }
}

void CartGrid::fill_comm_cell_lists(std::vector<int> &v,
                                    const Vec3i &lc,
                                    const Vec3i &hc)
{
    Vec3i c;
    for (c[0] = lc[0]; c[0] <= hc[0]; c[0]++) {
        for (c[1] = lc[1]; c[1] <= hc[1]; c[1]++) {
            for (c[2] = lc[2]; c[2] <= hc[2]; c[2]++) {
                v.push_back(linearize(c));
            }
        }
    }
}

void CartGrid::prepare_communication()
{
    m_exdescs.clear();
    m_exdescs.resize(n_neighbors());

    // Exclude first offset (0, 0, 0).
    for (const auto &offset : boost::make_iterator_range(
             std::next(std::begin(util::NeighborOffsets3D::raw)),
             std::end(util::NeighborOffsets3D::raw))) {
        Vec3i lc, hc;
        auto rank = proc_offset_to_rank(offset);

        // Skip self communication if it is unnecessary
        if (!self_comm_necessary() && rank == comm.rank())
            continue;

        // Send
        {
            auto i = neighbor_idx(rank);
            m_exdescs[i].dest = rank;
            std::tie(lc, hc)
                = impl::determine_send_receive_bounds(offset, 0, m_grid_size);
            fill_comm_cell_lists(m_exdescs[i].send, lc, hc);
        }

        // Receive in opposite direction. Otherwise send and receive order on
        // the processes won't match.
        {
            Vec3i opposite = {{-offset[0], -offset[1], -offset[2]}};
            auto neigh = proc_offset_to_rank(opposite);
            auto i = neighbor_idx(neigh);
            m_exdescs[i].dest = neigh;
            std::tie(lc, hc)
                = impl::determine_send_receive_bounds(opposite, 1, m_grid_size);
            fill_comm_cell_lists(m_exdescs[i].recv, lc, hc);
        }
    }
}

CartGrid::CartGrid(const boost::mpi::communicator &comm,
                   Vec3d box_size,
                   double min_cell_size)
    : ParallelLCGrid(comm, box_size, min_cell_size)
{
    create_grid();
    create_index_permutations();
    fill_neighranks();
    prepare_communication();
}

lidx CartGrid::n_local_cells()
{
    return m_grid_size[0] * m_grid_size[1] * m_grid_size[2];
}

gidx CartGrid::n_ghost_cells()
{
    int ggs
        = m_ghost_grid_size[0] * m_ghost_grid_size[1] * m_ghost_grid_size[2];
    return ggs - n_local_cells();
}

nidx CartGrid::n_neighbors()
{
    return m_neighranks.size();
}

rank CartGrid::neighbor_rank(nidx i)
{
    return m_neighranks[i];
}

lgidx CartGrid::cell_neighbor_index(lidx cellidx, int neigh)
{
    auto c = unlinearize(cellidx);
    auto nc = util::vadd_mod(c, util::NeighborOffsets3D::raw[neigh],
                             m_ghost_grid_size);
    return linearize(nc);
}

lgidx CartGrid::linearize(Vec3i c)
{
    auto idx = util::linearize(c, m_ghost_grid_size);
    return m_to_pargrid_order[idx];
}

Vec3i CartGrid::unlinearize(lgidx cidx)
{
    auto idx = m_from_pargrid_order[cidx];
    return util::unlinearize(idx, m_ghost_grid_size);
}

std::vector<GhostExchangeDesc> CartGrid::get_boundary_info()
{
    return m_exdescs;
}

lidx CartGrid::position_to_cell_index(Vec3d pos)
{
    Vec3i c;

    for (int i = 0; i < 3; ++i) {
        // Transform to process local coordinates
        double tpos = pos[i] - m_lowerleft[i];
        if (tpos < 0.0 || tpos >= m_localbox[i])
            throw std::domain_error("Particle not in local box");

        c[i] = tpos * m_inv_cell_size[i] + 1; // +1 to skip the ghost cells
    }
    return linearize(c);
}

rank CartGrid::position_to_rank(Vec3d pos)
{
    Vec3i proc;
    for (int i = 0; i < 3; ++i)
        proc[i]
            = static_cast<int>(pos[i] * m_inv_cell_size[i]) / m_grid_size[i];

    int rank;
    MPI_Cart_rank(comm_cart, proc.data(), &rank);
    return rank;
}

nidx CartGrid::neighbor_idx(rank r)
{
    // Search this rank in the local neighbor list and return its index
    // Use std::find here as 1) m_neighranks might not be sorted and 2) it has
    // at most 26 entries, so sequential search might not hurt that much.
    auto it = std::find(std::begin(m_neighranks), std::end(m_neighranks), r);
    if (*it != r)
        throw std::domain_error("Rank not a neighbor.");

    return std::distance(std::begin(m_neighranks), it);
}

nidx CartGrid::position_to_neighidx(Vec3d pos)
{
    // Determine the neighbor rank for locally known cell
    // Using position_to_rank here as it is the simpler code. Could also
    // search the neighboring cells of the cell where pos lies in.
    auto rank = position_to_rank(pos);
    return neighbor_idx(rank);
}

Vec3d CartGrid::cell_size()
{
    return m_cell_size;
}

Vec3i CartGrid::grid_size()
{
    Vec3i gs;
    // "m_grid_size" is the local number of cells.
    // Each process, however, has the same number of local cells per dimension.
    for (size_t i = 0; i < gs.size(); ++i) {
        gs[i] = m_grid_size[i] * node_grid[i];
    }
    return gs;
}

bool CartGrid::repartition(CellMetric m,
                           CellCellMetric ccm,
                           Thunk exchange_start_callback)
{
    UNUSED(m);
    UNUSED(ccm);
    UNUSED(exchange_start_callback);
    return false;
}

int CartGrid::global_hash(lgidx cellidx)
{
    Vec3i idx3d = unlinearize(cellidx);
    Vec3i dom = grid_size();
    Vec3i offset;

    for (size_t i = 0; i < idx3d.size(); ++i)
        offset[i] = m_procgrid_pos[i] * m_grid_size[i] - 1; // -1 for ghosts

    auto glo3didx = util::vadd_mod(idx3d, offset, dom);
    return util::linearize(glo3didx, dom);
}

} // namespace grids
} // namespace repa
