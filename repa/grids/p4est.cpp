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

//#ifdef HAVE_P4EST

#include <numeric>

#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>

#include "p4est.hpp"

#include "_compat.hpp"
#include "util/neighbor_offsets.hpp"

#ifdef __BMI2__
#include <x86intrin.h>
#endif

namespace repa {
namespace grids {

namespace impl {

RepartState::RepartState(const boost::mpi::communicator &comm_cart)
    : comm_cart(comm_cart),
      after_repart(false),
      nquads_per_proc(comm_cart.size())
{
}

void RepartState::reset()
{
    after_repart = false;
    std::fill(std::begin(nquads_per_proc), std::end(nquads_per_proc),
              static_cast<p4est_locidx_t>(0));
}

void RepartState::inc_nquads(rank proc)
{
    nquads_per_proc[proc]++;
}

void RepartState::allreduce()
{
    MPI_Allreduce(MPI_IN_PLACE, nquads_per_proc.data(), comm_cart.size(),
                  P4EST_MPI_LOCIDX, MPI_SUM, comm_cart);
}

// Returns the number of trailing zeros in an integer x.
static inline int count_trailing_zeros(int x)
{
    int z = 0;
    for (; (x & 1) == 0; x >>= 1)
        z++;
    return z;
}

/** Returns a bitset representing the boundary information of cell "idx".
 * Info is stored this way:
 * Bit 0  / Bit 1   / Bit 2  / Bit 3   / Bit 4  / Bit 5
 * X left / X right / Y left / Y right / Z left / Z right
 * Ordering: LSB ... MSB
 */
static int local_boundary_bitset(Vec3i idx, const Vec3i &grid_size)
{
    int ret = 0;
    for (int i = 0; i < 3; ++i) {
        if (PERIODIC(i) && idx[i] == 0)
            ret |= 1 << (2 * i);
        if (PERIODIC(i) && idx[i] == grid_size[i] - 1)
            ret |= 1 << (2 * i + 1);
    }
    return ret;
}

// Returns a global SFC-curve index for a given cell.
// Note: This is a global index on the Z-curve and not a local cell index to
// cells.
static gidx cell_morton_idx(Vec3i idx)
{
#ifdef __BMI2__
    static constexpr unsigned mask_x = 0x49249249;
    static constexpr unsigned mask_y = 0x92492492;
    static constexpr unsigned mask_z = 0x24924924;
    return _pdep_u32(idx[0], mask_x) | _pdep_u32(idx[1], mask_y)
           | _pdep_u32(idx[2], mask_z);
#else
    gidx res = 0;
    int res_bit = 1;

    for (int bit = 0; bit < 21; ++bit) {
        int mask = 1 << bit;
        for (int i = 0; i < 3; ++i) {
            if (idx[i] & mask)
                res |= res_bit;
            res_bit <<= 1;
        }
    }
    return res;
#endif
}

// Maps a position to the cartesian grid and returns the morton index of this
// coordinates
// Note: This is a global index on the Z-curve and not a local cell index to
// cells.
static gidx pos_morton_idx(Vec3d pos, const Vec3d &inv_cell_size)
{
    const Vec3i idx = {static_cast<int>(pos[0] * inv_cell_size[0]),
                       static_cast<int>(pos[1] * inv_cell_size[1]),
                       static_cast<int>(pos[2] * inv_cell_size[2])};
    return impl::cell_morton_idx(idx);
}

static inline Vec3i coord_to_cellindex(Vec3d coord, int tree_level)
{
    int scaling = 1 << tree_level;
    Vec3i idx;
    for (size_t i = 0; i < 3; ++i)
        idx[i] = coord[i] * scaling;
    return idx;
}

/** Scale coordinate coord, where 0.0 <= coord[i] < 1.0
 * to a cell index where 0 <= result[i] < maxidx,
 * where maxidx is determined by the tree/grid level.
 */
static inline Vec3i coord_to_cellindex(Vec3d coord, const Vec3i &tree_level)
{
    Vec3i idx;
    for (size_t i = 0; i < 3; ++i)
        idx[i] = coord[i] * (1 << tree_level[i]);
    return idx;
}

} // namespace impl

// Compute the grid- and bricksize according to box_l and maxrange
void P4estGrid::set_optimal_cellsize()
{
    // Compute number of cells and the cell size
    for (size_t i = 0; i < 3; ++i) {
        if (max_range > ROUND_ERROR_PREC * box_l[i])
            m_grid_size[i] = std::max<int>(box_l[i] / max_range, 1);
        else
            m_grid_size[i] = 1;

        m_cell_size[i] = box_l[i] / m_grid_size[i];
        m_inv_cell_size[i] = 1.0 / m_cell_size[i];
    }

    // Set number of trees to biggest common power of 2 of all dimensions
    m_grid_level = impl::count_trailing_zeros(m_grid_size[0] | m_grid_size[1]
                                              | m_grid_size[2]);
    for (size_t i = 0; i < 3; ++i)
        m_brick_size[i] = m_grid_size[i] >> m_grid_level;
}

void P4estGrid::create_grid()
{
    set_optimal_cellsize();

    if (!m_repartstate.after_repart) {
        // Keep old connectivity as the p4est destructor needs it
        auto oldconn = std::move(m_p8est_conn);
        // Create p8est structures
        m_p8est_conn = std::unique_ptr<p8est_connectivity_t>(
            p8est_connectivity_new_brick(m_brick_size[0], m_brick_size[1],
                                         m_brick_size[2], PERIODIC(0),
                                         PERIODIC(1), PERIODIC(2)));
        m_p8est = std::unique_ptr<p8est_t>(
            p8est_new_ext(comm_cart, m_p8est_conn.get(), 0, m_grid_level, true,
                          0, NULL, NULL));
    }

    // Information about first quads of each node
    // Assemble this as early as possible as it is necessary for
    // position_to_rank. As soon as this information is ready, we can start
    // migrating particles.
    m_node_first_cell_idx.resize(comm_cart.size() + 1);
    for (int i = 0; i < comm_cart.size(); ++i) {
        p8est_quadrant_t *q = &m_p8est->global_first_position[i];
        Vec3d xyz;
        p8est_qcoord_to_vertex(m_p8est_conn.get(), q->p.which_tree, q->x, q->y,
                               q->z, xyz.data());
        m_node_first_cell_idx[i] = impl::cell_morton_idx(
            impl::coord_to_cellindex(xyz, m_grid_level));
    }

    // Total number of quads
    int tmp = 1 << m_grid_level;
    while (tmp < m_grid_size[0] || tmp < m_grid_size[1] || tmp < m_grid_size[2])
        tmp <<= 1;
    m_node_first_cell_idx[comm_cart.size()] = tmp * tmp * tmp;

    if (m_repartstate.after_repart)
        m_repartstate.exchange_start_callback();

    auto p8est_ghost = std::unique_ptr<p8est_ghost_t>(
        p8est_ghost_new(m_p8est.get(), P8EST_CONNECT_CORNER));
    auto p8est_mesh = std::unique_ptr<p8est_mesh_t>(p8est_mesh_new_ext(
        m_p8est.get(), p8est_ghost.get(), 1, 1, 0, P8EST_CONNECT_CORNER));

    m_num_local_cells = m_p8est->local_num_quadrants;
    m_num_ghost_cells = p8est_ghost->ghosts.elem_count;
    int num_cells = m_num_local_cells + m_num_ghost_cells;

    std::unique_ptr<sc_array_t> ni
        = std::unique_ptr<sc_array_t>(sc_array_new(sizeof(int)));
    // Collect info about local cells
    m_p8est_shell.clear(); // Need to clear because we push_back
    m_global_idx.clear();
    m_p8est_shell.reserve(num_cells);
    for (int i = 0; i < m_num_local_cells; ++i) {
        p8est_quadrant_t *q
            = p8est_mesh_get_quadrant(m_p8est.get(), p8est_mesh.get(), i);
        Vec3d xyz;
        p8est_qcoord_to_vertex(m_p8est_conn.get(), p8est_mesh->quad_to_tree[i],
                               q->x, q->y, q->z, xyz.data());

        Vec3i idx = impl::coord_to_cellindex(
            xyz,
            p8est_tree_array_index(m_p8est->trees, p8est_mesh->quad_to_tree[i])
                ->maxlevel);

        // Cell on domain boundaries?
        int bndry = impl::local_boundary_bitset(idx, m_grid_size);
        m_p8est_shell.emplace_back(i, comm_cart.rank(),
                                   bndry ? impl::CellType::boundary
                                         : impl::CellType::inner,
                                   bndry, idx[0], idx[1], idx[2]);
        m_global_idx.push_back(
            impl::cell_morton_idx(impl::coord_to_cellindex(xyz, m_grid_level)));

        // Neighborhood
        for (int n = 0; n < 26; ++n) {
            m_p8est_shell[i].neighbor[n] = -1;
            p8est_mesh_get_neighbors(m_p8est.get(), p8est_ghost.get(),
                                     p8est_mesh.get(), i, n, NULL, NULL,
                                     ni.get());
            // Fully periodic, regular grid.
            if (ni->elem_count != 1)
                throw std::runtime_error("Error in neighborhood search.");

            int neighrank = *(int *)sc_array_index_int(ni.get(), 0);
            m_p8est_shell[i].neighbor[n] = neighrank;

            if (neighrank >= m_p8est->local_num_quadrants) {
                // Ghost cell on inner subdomain boundaries
                m_p8est_shell[i].shell = impl::CellType::boundary;
            }
            sc_array_truncate(ni.get());
        }
    }

    // Collect info about ghost cells
    for (int g = 0; g < m_num_ghost_cells; ++g) {
        p8est_quadrant_t *q
            = p8est_quadrant_array_index(&p8est_ghost->ghosts, g);
        Vec3d xyz;
        p8est_qcoord_to_vertex(m_p8est_conn.get(), q->p.piggy3.which_tree, q->x,
                               q->y, q->z, xyz.data());

        Vec3i idx = impl::coord_to_cellindex(
            xyz, p8est_tree_array_index(m_p8est->trees, q->p.piggy3.which_tree)
                     ->maxlevel);

        m_p8est_shell.emplace_back(g, p8est_mesh->ghost_to_proc[g],
                                   impl::CellType::ghost, 0, idx[0], idx[1],
                                   idx[2]);
        m_global_idx.push_back(
            impl::cell_morton_idx(impl::coord_to_cellindex(xyz, m_grid_level)));
    }
}

void P4estGrid::prepare_communication()
{
    int num_cells = n_local_cells() + n_ghost_cells();
    // List of cell indices for each process for send/recv
    std::vector<std::vector<lidx>> send_idx(comm_cart.size());
    std::vector<std::vector<gidx>> recv_idx(comm_cart.size());

    // Find all cells to be sent or received
    for (int i = 0; i < num_cells; ++i) {
        // Ghost cell? -> add to recv lists
        if (m_p8est_shell[i].shell == impl::CellType::ghost) {
            int nrank = m_p8est_shell[i].rank;
            if (nrank >= 0)
                recv_idx[nrank].push_back(i);
        }
        // Boundary cell? -> add to send lists
        if (m_p8est_shell[i].shell == impl::CellType::boundary) {
            // Add to all possible neighbors
            for (int n = 0; n < 26; ++n) {
                int nidx = m_p8est_shell[i].neighbor[n];
                // Invalid neighbor?
                if (nidx < 0
                    || m_p8est_shell[nidx].shell != impl::CellType::ghost)
                    continue;

                int nrank = m_p8est_shell[nidx].rank;
                if (nrank < 0)
                    continue;

                // Several neighbors n can be on the same process, therefore we
                // have to pay attention to add it only once.
                if (send_idx[nrank].empty() || send_idx[nrank].back() != i)
                    send_idx[nrank].push_back(i);
            }
        }
    }

    // Count the number of communications and assign indices for the
    // communication
    int num_comm_proc = 0;
    std::vector<int> comm_proc(comm_cart.size(), -1);
    for (int i = 0; i < comm_cart.size(); ++i) {
        if (send_idx[i].size() != 0 && recv_idx[i].size() != 0)
            comm_proc[i] = num_comm_proc++;
        else if (!(send_idx[i].size() == 0 && recv_idx[i].size() == 0))
            throw std::runtime_error(
                "Unexpected mismatch in send and receive lists.\n");
    }

    // Assemble ghost exchange descriptors
    m_exdescs.resize(num_comm_proc);
    m_neighranks.resize(num_comm_proc);
    for (int n = 0; n < comm_cart.size(); ++n) {
        if (comm_proc[n] == -1)
            continue;
        int index = comm_proc[n];
        m_neighranks[index] = n;
        m_exdescs[index].dest = n;
        m_exdescs[index].recv = std::move(recv_idx[n]);
        m_exdescs[index].send = std::move(send_idx[n]);
    }
}

void P4estGrid::reinitialize()
{
    create_grid();
    prepare_communication();
}

P4estGrid::P4estGrid(const boost::mpi::communicator &comm,
                     Vec3d box_size,
                     double min_cell_size)
    : ParallelLCGrid(comm, box_size, min_cell_size), m_repartstate(comm_cart)
{
    reinitialize();
}

lidx P4estGrid::n_local_cells()
{
    return m_num_local_cells;
}

gidx P4estGrid::n_ghost_cells()
{
    return m_num_ghost_cells;
}

nidx P4estGrid::n_neighbors()
{
    return m_neighranks.size();
}

rank P4estGrid::neighbor_rank(nidx i)
{
    if (i < 0 || i > m_neighranks.size())
        throw std::domain_error("Neighbor rank out of bounds.");
    return m_neighranks[i];
}

lgidx P4estGrid::cell_neighbor_index(lidx cellidx, fs_neighidx neigh)
{
    // Indices of the half shell neighbors in m_p8est_shell
    static const std::array<int, 27> to_p4est_order
    // Half-shell neighbors
        = {{-1, 1, 16, 3, 17, 22, 8, 23, 12, 5, 13, 24, 9, 25,
    //      ^^ unused. Cell itself is not stored in m_p8est_shell.
    // Rest (Full-shell \ Half-shell)
            0, 2, 4, 6, 7, 10, 11, 14, 15, 18, 19, 20, 21}};

    if (cellidx < 0 || cellidx >= n_local_cells())
        throw std::domain_error("Cell index outside of local subdomain");

    if (util::neighbor_type(neigh) == util::NeighborCellType::SELF)
        return cellidx;
    else
        return m_p8est_shell[cellidx].neighbor[to_p4est_order[neigh]];
}

std::vector<GhostExchangeDesc> P4estGrid::get_boundary_info()
{
    return m_exdescs;
}

lidx P4estGrid::position_to_cell_index(Vec3d pos)
{
    auto shellidxcomp = [](const impl::LocalShell &s, int idx) {
        int64_t sidx = impl::cell_morton_idx(s.coord);
        return sidx < idx;
    };

    auto needle
        = impl::pos_morton_idx(pos, m_inv_cell_size);

    auto shell_local_end = std::begin(m_p8est_shell) + n_local_cells();
    auto it
        = std::lower_bound(std::begin(m_p8est_shell),
                           // Only take into account local cells!
                           // This cannot be extended to ghost cells as these
                           // are not stored in SFC order in m_p8est_shell.
                           shell_local_end, needle, shellidxcomp);

    if (it != shell_local_end &&
        // Exclude finding cell 0 (lower_bound) if 0 is not the wanted result
        impl::cell_morton_idx(it->coord) == needle)
        return std::distance(std::begin(m_p8est_shell), it);
    else
        throw std::domain_error("Pos not in local domain.");
}

rank P4estGrid::position_to_rank(Vec3d pos)
{
    // Cell of pos might not be known on this process (not in m_p8est_shell).
    // Therefore, use the global first cell indices.
    auto it = std::upper_bound(
        std::begin(m_node_first_cell_idx), std::end(m_node_first_cell_idx),
        impl::pos_morton_idx(pos, m_inv_cell_size),
        [](int i, int64_t idx) { return i < idx; });

    return std::distance(std::begin(m_node_first_cell_idx), it) - 1;
}

nidx P4estGrid::position_to_neighidx(Vec3d pos)
{
    // Determine the neighbor rank for locally known cell
    // Using position_to_rank here as it is the simpler code. Could also
    // search the neighboring cells of the cell where pos lies in.
    auto rank = position_to_rank(pos);

    // Search this rank in the local neighbor list and return its index
    auto it = std::lower_bound(std::begin(m_neighranks), std::end(m_neighranks),
                               rank);
    if (it == std::end(m_neighranks) || *it != rank)
        throw std::domain_error("Position not in a ghost cell.");
    return std::distance(std::begin(m_neighranks), it);
}

Vec3d P4estGrid::cell_size()
{
    return m_cell_size;
}

Vec3i P4estGrid::grid_size()
{
    return m_grid_size;
}

bool P4estGrid::repartition(CellMetric m,
                            CellCellMetric ccm,
                            Thunk exchange_start_callback)
{
    UNUSED(ccm);
    // If this method exits early, successive calls to reinitialize() will
    // partition the grid uniformly.
    m_repartstate.reset();

    std::vector<double> weights = m();
    if (weights.size() != n_local_cells()) {
        throw std::runtime_error(
            "Metric only supplied " + std::to_string(weights.size())
            + "weights. Necessary: " + std::to_string(n_local_cells()));
    }

    // Determine prefix and target load
    double localsum
        = std::accumulate(std::begin(weights), std::end(weights), 0.0);
    double sum, prefix = 0; // Initialization is necessary on rank 0!
    MPI_Allreduce(&localsum, &sum, 1, MPI_DOUBLE, MPI_SUM, comm_cart);
    MPI_Exscan(&localsum, &prefix, 1, MPI_DOUBLE, MPI_SUM, comm_cart);
    double target = sum / comm_cart.size();

    // Determine new process boundaries in local subdomain
    // Evaluated for its side effect of setting part_nquads.
    std::accumulate(std::begin(weights), std::end(weights), prefix,
                    [this, target](double cellpref, double weight) {
                        int proc = std::min<int>(cellpref / target,
                                                 comm_cart.size() - 1);
                        m_repartstate.inc_nquads(proc);
                        return cellpref + weight;
                    });

    m_repartstate.allreduce();

    // TODO: Could try to steal quads from neighbors.
    //       Global reshifting (i.e. stealing from someone else than the direct
    //       neighbors) is not a good idea since it globally changes the metric.
    //       Anyways, this is most likely due to a bad quad/proc quotient.
    if (m_repartstate.nquads_per_proc[comm_cart.rank()] == 0) {
        fprintf(
            stderr,
            "[%i] No quads assigned to me. Cannot guarantee to work. Exiting\n",
            comm_cart.rank());
        fprintf(stderr,
                "[%i] Try changing the metric or reducing the number of "
                "processes\n",
                comm_cart.rank());
        errexit();
    }

    // Reinitialize the grid and prepare its internal datastructures for
    // querying by generic_dd.
    m_repartstate.after_repart = true;
    m_repartstate.exchange_start_callback = exchange_start_callback;

    p8est_partition_given(m_p8est.get(), m_repartstate.nquads_per_proc.data());
    reinitialize();

    return true;
}

int P4estGrid::global_hash(lgidx cellidx)
{
    return m_global_idx.at(cellidx);
}

} // namespace grids
} // namespace repa

//#endif // HAVE_P4EST
