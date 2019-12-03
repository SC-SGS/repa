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
#include "util/mpi_cart.hpp"
#include "util/neighbor_offsets.hpp"
#include "util/push_back_unique.hpp"
#include "util/vadd.hpp"
#include "util/vec_arith.hpp"

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
    using namespace util::vector_arithmetic;
    return any(c == 0) || any(c.as_expr() == m_ghost_grid_size - 1);
}

rank_type CartGrid::proc_offset_to_rank(const Vec3i &offset)
{
    using namespace util::vector_arithmetic;
    const Vec3i neighpos = (m_procgrid_pos + offset) % m_procgrid;
    return util::mpi_cart_rank(comm_cart, neighpos);
}

bool CartGrid::self_comm_necessary()
{
    using namespace util::vector_arithmetic;
    return any(m_procgrid < 2);
}

void CartGrid::fill_neighranks()
{
    m_neighranks.clear();

    for (const auto &offset : util::NeighborOffsets3D::raw) {
        // Push back unique neighbor ranks into m_neighbors
        const auto rank = proc_offset_to_rank(offset);
        if (rank != comm.rank() || self_comm_necessary())
            util::push_back_unique(m_neighranks, rank);
    }
}

void CartGrid::create_index_permutations()
{
    local_or_ghost_cell_index_type ncells = n_local_cells() + n_ghost_cells();
    m_to_pargrid_order.resize(ncells);
    m_from_pargrid_order.resize(ncells);

    local_cell_index_type localidx = 0;
    local_or_ghost_cell_index_type ghostidx = n_local_cells();
    for (local_or_ghost_cell_index_type i = 0; i < ncells; ++i) {
        const auto c = util::unlinearize(i, m_ghost_grid_size);
        local_or_ghost_cell_index_type &idx
            = is_ghost_cell(c) ? ghostidx : localidx;
        m_from_pargrid_order[idx] = i;
        m_to_pargrid_order[i] = idx++;
    }
}

void CartGrid::create_grid()
{
    using namespace util::vector_arithmetic;
    // Copy infos from deprecated pargrid variables
    m_procgrid = node_grid;
    m_procgrid_pos = node_pos;

    // Local box info
    m_localbox = box_l / m_procgrid;
    m_lowerleft = m_localbox * node_pos;

    if (max_range > ROUND_ERROR_PREC * box_l[0])
        m_grid_size = static_cast_vec<Vec3i>(m_localbox / max_range);
    else
        m_grid_size = constant_vec3(1);
    m_ghost_grid_size = m_grid_size + 2;

    m_cell_size = m_localbox / m_grid_size;
    m_inv_cell_size = 1.0 / m_cell_size;
}

void CartGrid::fill_comm_cell_lists(
    std::vector<local_or_ghost_cell_index_type> &v,
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
        const auto rank = proc_offset_to_rank(offset);

        // Skip self communication if it is unnecessary
        if (rank == comm_cart.rank() && !self_comm_necessary())
            continue;

        // Send
        {
            const auto i = neighbor_idx(rank);
            m_exdescs[i].dest = rank;
            std::tie(lc, hc)
                = impl::determine_send_receive_bounds(offset, 0, m_grid_size);
            fill_comm_cell_lists(m_exdescs[i].send, lc, hc);
        }

        // Receive in opposite direction. Otherwise send and receive order on
        // the processes won't match.
        {
            Vec3i opposite{-offset[0], -offset[1], -offset[2]};
            const auto neigh = proc_offset_to_rank(opposite);
            const auto i = neighbor_idx(neigh);
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

local_cell_index_type CartGrid::n_local_cells()
{
    return m_grid_size[0] * m_grid_size[1] * m_grid_size[2];
}

ghost_cell_index_type CartGrid::n_ghost_cells()
{
    const int ggs
        = m_ghost_grid_size[0] * m_ghost_grid_size[1] * m_ghost_grid_size[2];
    return ggs - n_local_cells();
}

rank_index_type CartGrid::n_neighbors()
{
    return m_neighranks.size();
}

rank_type CartGrid::neighbor_rank(rank_index_type i)
{
    return m_neighranks[i];
}

local_or_ghost_cell_index_type
CartGrid::cell_neighbor_index(local_cell_index_type cellidx, fs_neighidx neigh)
{
    using namespace util::vector_arithmetic;
    const auto coord = unlinearize(cellidx);
    const Vec3i neighbor_coord
        = (coord + util::NeighborOffsets3D::raw[neigh]) % m_ghost_grid_size;
    return linearize(neighbor_coord);
}

local_or_ghost_cell_index_type CartGrid::linearize(Vec3i c)
{
    const auto idx = util::linearize(c, m_ghost_grid_size);
    return m_to_pargrid_order[idx];
}

Vec3i CartGrid::unlinearize(local_or_ghost_cell_index_type cidx)
{
    const auto idx = m_from_pargrid_order[cidx];
    return util::unlinearize(idx, m_ghost_grid_size);
}

std::vector<GhostExchangeDesc> CartGrid::get_boundary_info()
{
    return m_exdescs;
}

local_cell_index_type CartGrid::position_to_cell_index(Vec3d pos)
{
    using namespace util::vector_arithmetic;
    const Vec3d relative_pos = pos - m_lowerleft;

    if (any(relative_pos < 0.0)
        || any(relative_pos.as_expr() >= m_localbox.as_expr()))
        throw std::domain_error("Particle not in local box");

    return linearize(static_cast_vec<Vec3i>(relative_pos * m_inv_cell_size)
                     + 1); // +1 to skip the ghost cells
}

rank_type CartGrid::position_to_rank(Vec3d pos)
{
    using namespace util::vector_arithmetic;
    return util::mpi_cart_rank(
        comm_cart, static_cast_vec<Vec3i>(pos * m_inv_cell_size) / m_grid_size);
}

rank_index_type CartGrid::neighbor_idx(rank_type r)
{
    // Search this rank in the local neighbor list and return its index
    // Use std::find here as 1) m_neighranks might not be sorted and 2) it has
    // at most 26 entries, so sequential search might not hurt that much.
    const auto it
        = std::find(std::begin(m_neighranks), std::end(m_neighranks), r);
    if (*it != r)
        throw std::domain_error("Rank not a neighbor.");

    return std::distance(std::begin(m_neighranks), it);
}

rank_index_type CartGrid::position_to_neighidx(Vec3d pos)
{
    // Determine the neighbor rank for locally known cell
    // Using position_to_rank here as it is the simpler code. Could also
    // search the neighboring cells of the cell where pos lies in.
    const auto rank = position_to_rank(pos);
    return neighbor_idx(rank);
}

Vec3d CartGrid::cell_size()
{
    return m_cell_size;
}

Vec3i CartGrid::grid_size()
{
    // "m_grid_size" is the local number of cells.
    // Each process, however, has the same number of local cells per dimension.
    using namespace util::vector_arithmetic;
    return m_grid_size * node_grid;
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

global_cell_index_type
CartGrid::global_hash(local_or_ghost_cell_index_type cellidx)
{
    // No need to define this away. Does currently not require extra data.
    using namespace util::vector_arithmetic;

    const Vec3i idx3d = unlinearize(cellidx);
    const Vec3i whole_domain = grid_size();
    const Vec3i prefix = m_procgrid_pos * m_grid_size - 1; // -1 for ghosts
    const Vec3i glo3didx = (idx3d + prefix) % whole_domain;

    return util::linearize(glo3didx, whole_domain);
}

} // namespace grids
} // namespace repa
