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

#include <numeric>

#include <functional>
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>

#include "p4est.hpp"

#include "_compat.hpp"
#include "util/neighbor_offsets.hpp"
#include "util/vec_arith.hpp"

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

void RepartState::inc_nquads(rank_type proc)
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

/** Returns true if "idx" is a boundary cell coordinate.
 */
static bool is_boundary(const Vec3i &idx, const Vec3i &grid_size)
{
    for (int i = 0; i < 3; ++i) {
        if (PERIODIC(i) && (idx[i] == 0 || idx[i] == grid_size[i] - 1))
            return true;
    }
    return false;
}

// Returns a global SFC-curve index for a given cell.
// Note: This is a global index on the Z-curve and not a local cell index to
// cells.
static global_cell_index_type cell_morton_idx(Vec3i idx)
{
#ifdef __BMI2__
    static constexpr unsigned mask_x = 0x49249249;
    static constexpr unsigned mask_y = 0x92492492;
    static constexpr unsigned mask_z = 0x24924924;
    return _pdep_u32(idx[0], mask_x) | _pdep_u32(idx[1], mask_y)
           | _pdep_u32(idx[2], mask_z);
#else
    global_cell_index_type res = 0;
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
static global_cell_index_type pos_morton_idx(Vec3d pos,
                                             const Vec3d &inv_cell_size)
{
    using namespace util::vector_arithmetic;
    return impl::cell_morton_idx(static_cast_vec<Vec3i>(pos * inv_cell_size));
}

static inline Vec3i coord_to_cellindex(Vec3d coord, int tree_level)
{
    using namespace util::vector_arithmetic;
    const int scaling = 1 << tree_level;
    return static_cast_vec<Vec3i>(coord * static_cast<double>(scaling));
}

/** Scale coordinate coord, where 0.0 <= coord[i] < 1.0
 * to a cell index where 0 <= result[i] < maxidx,
 * where maxidx is determined by the tree/grid level.
 */
// static inline Vec3i coord_to_cellindex(Vec3d coord, const Vec3i &tree_level)
//{
//    using namespace util::vector_arithmetic;
//    return static_cast_vec<Vec3i>(coord * (1 << tree_level));
//}

template <typename T>
T virtual_regular_grid_end(T grid_level, const Vec3<T> &grid_size)
{
    // Next largest power of 2
    using namespace repa::util::vector_arithmetic;
    T gs = 1 << grid_level;
    while (any(grid_size > gs))
        gs <<= 1;
    return gs * gs * gs;
}

static Vec3d quadrant_to_coords(const p8est_quadrant_t *q,
                                p8est_connectivity_t *p8est_conn)
{
    Vec3d xyz;
    p8est_qcoord_to_vertex(p8est_conn, q->p.which_tree, q->x, q->y, q->z,
                           xyz.data());
    return xyz;
}

static Vec3d quadrant_to_coords(const p8est_quadrant_t *q,
                                p8est_connectivity_t *p8est_conn,
                                p4est_topidx_t tree_id)
{
    Vec3d xyz;
    p8est_qcoord_to_vertex(p8est_conn, tree_id, q->x, q->y, q->z, xyz.data());
    return xyz;
}

} // namespace impl

// Compute the grid- and bricksize according to box_l and maxrange
void P4estGrid::set_optimal_cellsize()
{
    using namespace util::vector_arithmetic;
    // Compute number of cells and the cell size
    if (max_range > ROUND_ERROR_PREC * box_l[0])
        m_grid_size = static_cast_vec<Vec3i>(box_l / max_range);
    else
        m_grid_size = constant_vec3(1);

    m_cell_size = box_l / m_grid_size;
    m_inv_cell_size = 1.0 / m_cell_size;

    // Set number of trees to biggest common power of 2 of all dimensions
    m_grid_level
        = impl::count_trailing_zeros(m_grid_size.foldl(std::bit_or<>{}));
    m_brick_size = m_grid_size >> m_grid_level;
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
    m_node_first_cell_idx.clear();
    m_node_first_cell_idx.reserve(comm_cart.size() + 1);
    for (int i = 0; i < comm_cart.size(); ++i) {
        m_node_first_cell_idx.push_back(
            impl::cell_morton_idx(impl::coord_to_cellindex(
                impl::quadrant_to_coords(&m_p8est->global_first_position[i],
                                         m_p8est_conn.get()),
                m_grid_level)));
    }
    // Total number of quads of the virtual regular grid
    m_node_first_cell_idx.push_back(
        impl::virtual_regular_grid_end(m_grid_level, m_grid_size));

    if (m_repartstate.after_repart)
        m_repartstate.exchange_start_callback();

    auto p8est_ghost = std::unique_ptr<p8est_ghost_t>(
        p8est_ghost_new(m_p8est.get(), P8EST_CONNECT_CORNER));
    auto p8est_mesh = std::unique_ptr<p8est_mesh_t>(p8est_mesh_new_ext(
        m_p8est.get(), p8est_ghost.get(), 1, 1, 0, P8EST_CONNECT_CORNER));

    m_num_local_cells = m_p8est->local_num_quadrants;
    m_num_ghost_cells = p8est_ghost->ghosts.elem_count;
    const local_or_ghost_cell_index_type num_cells
        = m_num_local_cells + m_num_ghost_cells;

    std::unique_ptr<sc_array_t> ni
        = std::unique_ptr<sc_array_t>(sc_array_new(sizeof(int)));
    // Collect info about local cells
#ifdef GLOBAL_HASH_NEEDED
    m_global_idx.clear();
#endif
    m_p8est_cell_info.clear(); // Need to clear because we push_back
    m_p8est_cell_info.reserve(num_cells);
    for (local_cell_index_type i = 0; i < m_num_local_cells; ++i) {
        const Vec3d xyz = impl::quadrant_to_coords(
            p8est_mesh_get_quadrant(m_p8est.get(), p8est_mesh.get(),
                                    static_cast<p4est_locidx_t>(i)),
            m_p8est_conn.get(), p8est_mesh->quad_to_tree[i]);

        const Vec3i idx = impl::coord_to_cellindex(
            xyz,
            p8est_tree_array_index(m_p8est->trees, p8est_mesh->quad_to_tree[i])
                ->maxlevel);

        m_p8est_cell_info.emplace_back(comm_cart.rank(),
                                       impl::is_boundary(idx, m_grid_size)
                                           ? impl::CellType::boundary
                                           : impl::CellType::inner,
                                       idx);
#ifdef GLOBAL_HASH_NEEDED
        m_global_idx.push_back(
            impl::cell_morton_idx(impl::coord_to_cellindex(xyz, m_grid_level)));
#endif

        // Neighborhood
        for (int n = 0; n < 26; ++n) {
            p8est_mesh_get_neighbors(m_p8est.get(), p8est_ghost.get(),
                                     p8est_mesh.get(), i, n, NULL, NULL,
                                     ni.get());
            // Fully periodic, regular grid.
            assert(ni->elem_count == 1);

            const int neighidx = *(int *)sc_array_index_int(ni.get(), 0);
            m_p8est_cell_info[i].neighbor[n] = neighidx;

            if (neighidx >= m_p8est->local_num_quadrants) {
                // Ghost cell on inner subdomain boundaries
                m_p8est_cell_info[i].cell_type = impl::CellType::boundary;
            }
            sc_array_truncate(ni.get());
        }
    }

    // Collect info about ghost cells
    for (ghost_cell_index_type g = 0; g < m_num_ghost_cells; ++g) {
        const p8est_quadrant_t *q
            = p8est_quadrant_array_index(&p8est_ghost->ghosts, g);
        const Vec3d xyz = impl::quadrant_to_coords(q, m_p8est_conn.get(),
                                                   q->p.piggy3.which_tree);
        const Vec3i idx = impl::coord_to_cellindex(
            xyz, p8est_tree_array_index(m_p8est->trees, q->p.piggy3.which_tree)
                     ->maxlevel);

        m_p8est_cell_info.emplace_back(p8est_mesh->ghost_to_proc[g],
                                       impl::CellType::ghost, idx);
#ifdef GLOBAL_HASH_NEEDED
        m_global_idx.push_back(
            impl::cell_morton_idx(impl::coord_to_cellindex(xyz, m_grid_level)));
#endif
    }
}

void P4estGrid::prepare_communication()
{
    const local_or_ghost_cell_index_type num_cells
        = n_local_cells() + n_ghost_cells();
    // List of cell indices for each process for send/recv
    std::vector<std::vector<local_cell_index_type>> send_idx(comm_cart.size());
    std::vector<std::vector<ghost_cell_index_type>> recv_idx(comm_cart.size());

    // Find all cells to be sent or received
    for (local_or_ghost_cell_index_type i = 0; i < num_cells; ++i) {
        const auto &cell_info = m_p8est_cell_info[i];
        // Ghost cell? -> add to recv lists
        if (cell_info.cell_type == impl::CellType::ghost) {
            const rank_type nrank = m_p8est_cell_info[i].owner_rank;
            assert(nrank >= 0 && nrank < comm.size());
            recv_idx[nrank].push_back(i);
        }
        // Boundary cell? -> add to send lists
        else if (cell_info.cell_type == impl::CellType::boundary) {
            // Add to all possible neighbors
            for (local_or_ghost_cell_index_type nidx : cell_info.neighbor) {
                assert(nidx >= 0 && nidx < num_cells);
                // Only ghost cells hold possible neighboring processes
                if (m_p8est_cell_info[nidx].cell_type != impl::CellType::ghost)
                    continue;

                const rank_type nrank = m_p8est_cell_info[nidx].owner_rank;
                assert(nrank >= 0 && nrank < comm.size());

                // Several neighbor cells nidx can be on the same process.
                // We must only add it once.
                if (send_idx[nrank].empty() || send_idx[nrank].back() != i)
                    send_idx[nrank].push_back(i);
            }
        }
    }

    // Assemble ghost exchange information and neighbors
    m_exdescs.clear();
    m_neighranks.clear();
    for (rank_type neighrank = 0; neighrank < comm_cart.size(); ++neighrank) {

        if (send_idx[neighrank].empty() || recv_idx[neighrank].empty()) {
            assert(send_idx[neighrank].empty() && recv_idx[neighrank].empty());
            continue;
        }

        m_neighranks.push_back(neighrank);
        m_exdescs.emplace_back(neighrank, std::move(recv_idx[neighrank]),
                               std::move(send_idx[neighrank]));
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

local_cell_index_type P4estGrid::n_local_cells()
{
    return m_num_local_cells;
}

ghost_cell_index_type P4estGrid::n_ghost_cells()
{
    return m_num_ghost_cells;
}

rank_index_type P4estGrid::n_neighbors()
{
    return m_neighranks.size();
}

rank_type P4estGrid::neighbor_rank(rank_index_type i)
{
    assert(i >= 0 && i < n_neighbors());
    return m_neighranks[i];
}

local_or_ghost_cell_index_type
P4estGrid::cell_neighbor_index(local_cell_index_type cellidx, fs_neighidx neigh)
{
    // Indices of the half shell neighbors in m_p8est_cell_info
    static const std::array<int, 27> to_p4est_order
        // Half-shell neighbors
        = {{-1, 1, 16, 3, 17, 22, 8, 23, 12, 5, 13, 24, 9, 25,
            //      ^^ unused. Cell itself is not stored in m_p8est_cell_info.
            // Rest (Full-shell \ Half-shell)
            0, 2, 4, 6, 7, 10, 11, 14, 15, 18, 19, 20, 21}};

    assert(cellidx >= 0 && cellidx < n_local_cells());

    if (util::neighbor_type(neigh) == util::NeighborCellType::SELF)
        return cellidx;
    else
        return m_p8est_cell_info[cellidx].neighbor[to_p4est_order[neigh]];
}

std::vector<GhostExchangeDesc> P4estGrid::get_boundary_info()
{
    return m_exdescs;
}

local_cell_index_type P4estGrid::position_to_cell_index(Vec3d pos)
{
    auto shellidxcomp
        = [](const impl::CellInfo &s, global_cell_index_type idx) {
              return impl::cell_morton_idx(s.coord) < idx;
          };

    auto needle = impl::pos_morton_idx(pos, m_inv_cell_size);

    auto shell_local_end = std::begin(m_p8est_cell_info) + n_local_cells();
    auto it
        = std::lower_bound(std::begin(m_p8est_cell_info),
                           // Only take into account local cells!
                           // This cannot be extended to ghost cells as these
                           // are not stored in SFC order in m_p8est_cell_info.
                           shell_local_end, needle, shellidxcomp);

    if (it != shell_local_end &&
        // Exclude finding cell 0 (lower_bound) if 0 is not the wanted result
        impl::cell_morton_idx(it->coord) == needle)
        return std::distance(std::begin(m_p8est_cell_info), it);
    else
        throw std::domain_error("Pos not in local domain.");
}

rank_type P4estGrid::position_to_rank(Vec3d pos)
{
    // Cell of pos might not be known on this process (not in
    // m_p8est_cell_info). Therefore, use the global first cell indices.
    auto it = std::upper_bound(
        std::begin(m_node_first_cell_idx), std::end(m_node_first_cell_idx),
        impl::pos_morton_idx(pos, m_inv_cell_size),
        [](global_cell_index_type i, global_cell_index_type idx) {
            return i < idx;
        });

    return std::distance(std::begin(m_node_first_cell_idx), it) - 1;
}

rank_index_type P4estGrid::position_to_neighidx(Vec3d pos)
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

    const auto weights = m();
    assert(weights.size() == n_local_cells());

    // Determine prefix and target load
    const double localsum
        = std::accumulate(std::begin(weights), std::end(weights), 0.0);
    double sum, prefix = 0; // Initialization is necessary on rank 0!
    MPI_Allreduce(&localsum, &sum, 1, MPI_DOUBLE, MPI_SUM, comm_cart);
    MPI_Exscan(&localsum, &prefix, 1, MPI_DOUBLE, MPI_SUM, comm_cart);
    const double target = sum / comm_cart.size();

    // Determine new process boundaries in local subdomain
    // Evaluated for its side effect of setting part_nquads.
    std::accumulate(std::begin(weights), std::end(weights), prefix,
                    [this, target](double cellpref, double weight) {
                        rank_type proc = std::min<rank_type>(
                            cellpref / target, comm_cart.size() - 1);
                        m_repartstate.inc_nquads(proc);
                        return cellpref + weight;
                    });

    m_repartstate.allreduce();

    // TODO: Could try to steal quads from neighbors.
    //       Global reshifting (i.e. stealing from someone else than the direct
    //       neighbors) is not a good idea since it globally changes the metric.
    //       Anyways, this is most likely due to a bad quad/proc quotient.
    assert(m_repartstate.nquads_per_proc[comm_cart.rank()] > 0);

    // Reinitialize the grid and prepare its internal datastructures for
    // querying by generic_dd.
    m_repartstate.after_repart = true;
    m_repartstate.exchange_start_callback = exchange_start_callback;

    p8est_partition_given(m_p8est.get(), m_repartstate.nquads_per_proc.data());
    reinitialize();

    return true;
}

global_cell_index_type
P4estGrid::global_hash(local_or_ghost_cell_index_type cellidx)
{
#ifdef GLOBAL_HASH_NEEDED
    return m_global_idx.at(cellidx);
#else
    return 0;
#endif
}

} // namespace grids
} // namespace repa
