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

#include <algorithm>
#include <functional>
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_extended.h>
#include <vector>

#include "p4est.hpp"

#include "_compat.hpp"
#include "util/neighbor_offsets.hpp"
#include "util/p4est_deleter.hpp"
#include "util/range.hpp"
#include "util/vec_arith.hpp"

// clang-format off
#ifdef __INTEL_COMPILER
#   ifdef __BMI__
#       include <immintrin.h>
#       define HAS_PDEP_INTRINSIC
#   endif
#else
#   ifdef __BMI2__
#       include <x86intrin.h>
#       define HAS_PDEP_INTRINSIC
#   endif
#endif
// clang-format on

namespace repa {
namespace grids {

namespace impl {

void init_p4est_logging()
{
    static bool initialized = false;
    if (!initialized) {
        p4est_init(NULL, SC_LP_ERROR);
        initialized = true;
    }
}

enum class CellType { inner = 0, boundary = 1, ghost = 2 };
struct CellInfo {
    const rank_type
        owner_rank;     // the rank of this cell (equals this_node for locals)
    CellType cell_type; // shell information (0: inner local cell, 1: boundary
                        // local cell, 2: ghost cell)
                        // boundary Cells with shell-type 2 are set if the are
                        // in the periodic halo
    std::array<local_or_ghost_cell_index_type, 26>
        neighbor; // unique index of the fullshell neighborhood cells (as in
                  // p8est); only 26 because cell itself is not included.
    const Vec3i coord; // cartesian coordinates of the cell
                       // For a globally unique index (on a virtual regular
                       // grid) call impl::cell_morton_index(coord).

    CellInfo(rank_type rank, CellType shell, const Vec3i &coord)
        : owner_rank(rank), cell_type(shell), coord(coord)
    {
    }
};

struct RepartState {
    const boost::mpi::communicator &comm_cart;
    bool after_repart;
    std::vector<p4est_locidx_t> nquads_per_proc;
    std::function<void()> exchange_start_callback;

    RepartState(const boost::mpi::communicator &comm_cart)
        : comm_cart(comm_cart),
          after_repart(false),
          nquads_per_proc(comm_cart.size())
    {
    }

    inline void reset()
    {
        after_repart = false;
        std::fill(std::begin(nquads_per_proc), std::end(nquads_per_proc),
                  static_cast<p4est_locidx_t>(0));
    }
    inline void inc_nquads(rank_type proc)
    {
        nquads_per_proc[proc]++;
    }
    inline void allreduce()
    {
        MPI_Allreduce(MPI_IN_PLACE, nquads_per_proc.data(), comm_cart.size(),
                      boost::mpi::get_mpi_datatype<decltype(
                          nquads_per_proc)::value_type>(),
                      MPI_SUM, comm_cart);
    }
};

// Returns the number of trailing zeros in an integer x.
static inline int count_trailing_zeros(int x)
{
    assert(x != 0);
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
#ifdef HAS_PDEP_INTRINSIC
    static constexpr unsigned mask_x = 0x49249249;
    static constexpr unsigned mask_y = 0x92492492;
    static constexpr unsigned mask_z = 0x24924924;
    return global_cell_index_type{_pdep_u32(idx[0], mask_x)
                                  | _pdep_u32(idx[1], mask_y)
                                  | _pdep_u32(idx[2], mask_z)};
#else
    global_cell_index_type::value_type res = 0;
    global_cell_index_type::value_type res_bit = 1;
    constexpr size_t max_nbits_per_dim = sizeof(global_cell_index_type) * 8 / 3;

    for (size_t bit = 0; bit < max_nbits_per_dim; ++bit) {
        int mask = 1 << bit;
        for (int i = 0; i < 3; ++i) {
            if (idx[i] & mask)
                res |= res_bit;
            res_bit <<= 1;
        }
    }
    return global_cell_index_type{res};
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

template <typename Ret, typename T>
Ret virtual_regular_grid_end(T grid_level, const Vec3<T> &grid_size)
{
    // Next largest power of 2
    using namespace repa::util::vector_arithmetic;
    Ret gs = 1 << grid_level;
    while (any(grid_size > gs))
        gs <<= 1;
    return Ret{gs * gs * gs};
}

static Vec3d quadrant_to_coords(const p8est_quadrant_t *q,
                                p8est_connectivity_t *p8est_conn)
{
    Vec3d xyz;
    assert(q->p.which_tree < p8est_conn->num_trees);
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

/** Fills "first_morton_index_per_proc" with indices of quadrants
 * along the SFC on a regular grid of size "grid_size" and level "grid_level".
 */
static void assemble_virtual_regular_grid_partitioning(
    std::vector<global_cell_index_type> &first_morton_index_per_proc,
    int grid_level,
    repa::Vec3i grid_size,
    rank_type nproc,
    p8est_t *p8est,
    p8est_connectivity_t *p8est_conn)
{
    first_morton_index_per_proc.clear();
    first_morton_index_per_proc.reserve(nproc + 1);
    for (int i = 0; i < nproc; ++i) {
        first_morton_index_per_proc.push_back(
            impl::cell_morton_idx(impl::coord_to_cellindex(
                impl::quadrant_to_coords(&p8est->global_first_position[i],
                                         p8est_conn),
                grid_level)));
    }
    // Total number of quads of the virtual regular grid
    first_morton_index_per_proc.push_back(global_cell_index_type{
        impl::virtual_regular_grid_end<global_cell_index_type::value_type>(
            grid_level, grid_size)});
}

} // namespace impl

struct P4estGrid::_P4estGrid_impl {
    int m_grid_level;

    // Number of grid cells in total and per tree.
    Vec3i m_grid_size, m_brick_size;
    // Cell size (box_size / m_grid_size)
    Vec3d m_cell_size, m_inv_cell_size;

    // p4est data structures
    std::unique_ptr<p8est_connectivity_t> m_p8est_conn;
    std::unique_ptr<p8est_t> m_p8est;
    local_cell_index_type m_num_local_cells;
    ghost_cell_index_type m_num_ghost_cells;

    // helper data structures
    std::vector<global_cell_index_type>
        m_node_first_cell_idx; // First morton index on a virtual regular grid
                               // spanning all processes.
#ifdef GLOBAL_HASH_NEEDED
    std::vector<global_cell_index_type>
        m_global_idx; //< Global virtual morton index of cells (for
                      // global_hash()) could be removed in non-test builds.
#endif
    std::vector<impl::CellInfo> m_p8est_cell_info; //< Raw cell info from p4est

    // comm data structures
    std::vector<GhostExchangeDesc> m_exdescs;
    std::vector<rank_type> m_neighranks;

    void set_optimal_cellsize();
    void create_grid();
    void init_grid_cells(p8est_ghost_t *p8est_ghost, p8est_mesh_t *p8est_mesh);
    void prepare_communication();
    // Reinitialized the grid (instantiation or after repartitioning)
    void reinitialize();
    bool repartition(CellMetric m,
                     CellCellMetric ccm,
                     Thunk exchange_start_callback);

    local_or_ghost_cell_index_type
    p4est_index_to_locghost(p4est_locidx_t idx) const
    {
        // P4est coding as follows:
        // 0 <= idx < n_local_cells ==> local cell index
        // n_local_cells <= idx < n_local_cells + n_ghost_cells ==> ghost
        // cell index
        assert(idx >= 0 && idx < m_num_local_cells + m_num_ghost_cells);

        if (idx >= m_num_local_cells)
            return ghost_cell_index_type{idx - m_num_local_cells};
        else
            return local_cell_index_type{idx};
    }

    p4est_locidx_t
    locghost_to_p4est_index(local_or_ghost_cell_index_type idx) const
    {
        p4est_locidx_t result;
        idx.visit(
            [&result
#ifndef NDEBUG
             ,
             this
#endif
        ](local_cell_index_type lidx) {
                result = static_cast<int>(lidx);
                assert(result >= 0 && result < m_num_local_cells);
            },
            [&result, this](ghost_cell_index_type gidx) {
                result = static_cast<int>(gidx) + m_num_local_cells;
                assert(result >= m_num_local_cells
                       && result < m_num_local_cells + m_num_ghost_cells);
            });
        return result;
    }

    impl::RepartState m_repartstate;

    const boost::mpi::communicator &comm_cart;
    const Vec3d box_size;
    const double min_cell_size;

    _P4estGrid_impl(const boost::mpi::communicator &comm_cart,
                    const Vec3d &box_size,
                    double min_cell_size)
        : m_repartstate(comm_cart),
          comm_cart(comm_cart),
          box_size(box_size),
          min_cell_size(min_cell_size)
    {
        impl::init_p4est_logging();
    }
};

// Compute the grid- and bricksize according to box_size and maxrange
void P4estGrid::_P4estGrid_impl::set_optimal_cellsize()
{
    using namespace util::vector_arithmetic;
    // Compute number of cells and the cell size
    m_grid_size = static_cast_vec<Vec3i>(box_size / min_cell_size);

    m_cell_size = box_size / m_grid_size;
    m_inv_cell_size = 1.0 / m_cell_size;

    // Set number of trees to biggest common power of 2 of all dimensions
    using bit_or = std::bit_or<decltype(
        m_grid_size)::value_type>; // For compatibility with old
                                   // standard libraries, that do not
                                   // support std::bit_or<T = void>
    m_grid_level = impl::count_trailing_zeros(m_grid_size.foldl(bit_or()));
    m_brick_size = m_grid_size >> m_grid_level;
}

void P4estGrid::_P4estGrid_impl::create_grid()
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
    impl::assemble_virtual_regular_grid_partitioning(
        m_node_first_cell_idx, m_grid_level, m_grid_size, comm_cart.size(),
        m_p8est.get(), m_p8est_conn.get());

    if (m_repartstate.after_repart)
        m_repartstate.exchange_start_callback();

    // Ghost and mesh are only necessary for init_grid_cells and can be
    // deleted afterwards.
    auto p8est_ghost = std::unique_ptr<p8est_ghost_t>(
        p8est_ghost_new(m_p8est.get(), P8EST_CONNECT_CORNER));
    auto p8est_mesh = std::unique_ptr<p8est_mesh_t>(p8est_mesh_new_ext(
        m_p8est.get(), p8est_ghost.get(), 1, 1, 0, P8EST_CONNECT_CORNER));

    m_num_local_cells = local_cell_index_type{m_p8est->local_num_quadrants};
    m_num_ghost_cells = ghost_cell_index_type{p8est_ghost->ghosts.elem_count};

    init_grid_cells(p8est_ghost.get(), p8est_mesh.get());
}

void P4estGrid::_P4estGrid_impl::init_grid_cells(p8est_ghost_t *p8est_ghost,
                                                 p8est_mesh_t *p8est_mesh)
{
    // "ni" is defined outside of all loops to avoid calling "sc_array_new" and
    // the corresponding destroy function on every loop iteration.
    std::unique_ptr<sc_array_t> ni
        = std::unique_ptr<sc_array_t>(sc_array_new(sizeof(int)));

#ifdef GLOBAL_HASH_NEEDED
    m_global_idx.clear();
#endif
    m_p8est_cell_info.clear();
    m_p8est_cell_info.reserve(m_num_local_cells + m_num_ghost_cells);
    for (const auto i : util::range(m_num_local_cells)) {
        const Vec3d xyz = impl::quadrant_to_coords(
            p8est_mesh_get_quadrant(m_p8est.get(), p8est_mesh,
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

        auto &cell_info = m_p8est_cell_info[i];
        // Collect neighborhood information
        for (int n = 0; n < 26; ++n) {
            p8est_mesh_get_neighbors(m_p8est.get(), p8est_ghost, p8est_mesh, i,
                                     n, NULL, NULL, ni.get());
            // Fully periodic, regular grid.
            assert(ni->elem_count == 1);

            cell_info.neighbor[n] = p4est_index_to_locghost(
                *(int *)sc_array_index_int(ni.get(), 0));

            sc_array_truncate(ni.get());
        }

        // Mark inner boundary cells
        if (std::any_of(std::begin(cell_info.neighbor),
                        std::end(cell_info.neighbor),
                        [](const local_or_ghost_cell_index_type &neighidx) {
                            return neighidx.is<ghost_cell_index_type>();
                        }))
            cell_info.cell_type = impl::CellType::boundary;
    }

    // Collect info about ghost cells
    for (const auto g : util::range(m_num_ghost_cells)) {
        const p8est_quadrant_t *q
            = p8est_quadrant_array_index(&p8est_ghost->ghosts, g);
        const Vec3d xyz = impl::quadrant_to_coords(q, m_p8est_conn.get(),
                                                   q->p.piggy3.which_tree);
        const Vec3i idx = impl::coord_to_cellindex(
            xyz, p8est_tree_array_index(m_p8est->trees, q->p.piggy3.which_tree)
                     ->maxlevel);

        assert(p8est_mesh->ghost_to_proc[g] >= 0
               && p8est_mesh->ghost_to_proc[g] < comm_cart.size());

        m_p8est_cell_info.emplace_back(p8est_mesh->ghost_to_proc[g],
                                       impl::CellType::ghost, idx);
#ifdef GLOBAL_HASH_NEEDED
        m_global_idx.push_back(
            impl::cell_morton_idx(impl::coord_to_cellindex(xyz, m_grid_level)));
#endif
    }
}

void P4estGrid::_P4estGrid_impl::prepare_communication()
{
    const p4est_locidx_t num_cells = m_num_local_cells + m_num_ghost_cells;

    // List of cell indices for each process for send/recv
    std::vector<std::vector<local_cell_index_type>> send_idx(comm_cart.size());
    std::vector<std::vector<ghost_cell_index_type>> recv_idx(comm_cart.size());

    // Find all cells to be sent or received
    for (p4est_locidx_t i = 0; i < num_cells; ++i) {
        const auto &cell_info = m_p8est_cell_info[i];

        switch (cell_info.cell_type) {
        case impl::CellType::ghost:
            // Collect receive information of ghost cells
            recv_idx[cell_info.owner_rank].push_back(
                p4est_index_to_locghost(i).as<ghost_cell_index_type>());
            break;
        case impl::CellType::boundary:
            // Collect send information of boundary cells (check all neighboring
            // cells for neighboring processes)
            for (const auto &neighcell : cell_info.neighbor) {
                const auto &neigh_info
                    = m_p8est_cell_info[locghost_to_p4est_index(neighcell)];
                if (neigh_info.cell_type != impl::CellType::ghost)
                    continue;

                // Several neighbor cells can be on the same process.
                // We must only add it once. If "i" several neighbor cells with
                // the same owner, "i" must have been added in the
                // current iteration of the outer loop, and, thus, last in the
                // vector.
                if (send_idx[neigh_info.owner_rank].empty()
                    || send_idx[neigh_info.owner_rank].back() != i)
                    send_idx[neigh_info.owner_rank].push_back(
                        p4est_index_to_locghost(i).as<local_cell_index_type>());
            }
            break;
        default:
            // Nothing. Only interested in boundary and ghost cells.
            break;
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

bool P4estGrid::_P4estGrid_impl::repartition(CellMetric m,
                                             CellCellMetric ccm,
                                             Thunk exchange_start_callback)
{
    UNUSED(ccm);
    // If this method exits early, successive calls to reinitialize() will
    // partition the grid uniformly.
    m_repartstate.reset();

    const auto weights = m();
    assert(weights.size() == static_cast<size_t>(m_num_local_cells));

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
    ensure(m_repartstate.nquads_per_proc[comm_cart.rank()] > 0,
           "Detected a process with no quads assigned.");

    // Reinitialize the grid and prepare its internal datastructures for
    // querying by generic_dd.
    m_repartstate.after_repart = true;
    m_repartstate.exchange_start_callback = exchange_start_callback;

    p8est_partition_given(m_p8est.get(), m_repartstate.nquads_per_proc.data());
    reinitialize();

    return true;
}

void P4estGrid::_P4estGrid_impl::reinitialize()
{
    create_grid();
    prepare_communication();
}

P4estGrid::P4estGrid(const boost::mpi::communicator &comm,
                     Vec3d box_size,
                     double min_cell_size)
    : ParallelLCGrid(comm, box_size, min_cell_size),
      _impl(
          std::make_unique<_P4estGrid_impl>(comm_cart, box_size, min_cell_size))

{
    _impl->reinitialize();
}

local_cell_index_type P4estGrid::n_local_cells() const
{
    return _impl->m_num_local_cells;
}

ghost_cell_index_type P4estGrid::n_ghost_cells() const
{
    return _impl->m_num_ghost_cells;
}

util::const_span<rank_type> P4estGrid::neighbor_ranks() const
{
    return util::make_const_span(_impl->m_neighranks);
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
        return _impl->m_p8est_cell_info[cellidx]
            .neighbor[to_p4est_order[neigh]];
}

util::const_span<GhostExchangeDesc> P4estGrid::get_boundary_info()
{
    return util::make_const_span(_impl->m_exdescs);
}

local_cell_index_type P4estGrid::position_to_cell_index(Vec3d pos)
{
    auto shellidxcomp
        = [](const impl::CellInfo &s, global_cell_index_type idx) {
              return impl::cell_morton_idx(s.coord) < idx;
          };

    auto needle = impl::pos_morton_idx(pos, _impl->m_inv_cell_size);

    auto shell_local_end
        = std::begin(_impl->m_p8est_cell_info) + n_local_cells();
    auto it
        = std::lower_bound(std::begin(_impl->m_p8est_cell_info),
                           // Only take into account local cells!
                           // This cannot be extended to ghost cells as these
                           // are not stored in SFC order in m_p8est_cell_info.
                           shell_local_end, needle, shellidxcomp);

    if (it != shell_local_end &&
        // Exclude finding cell 0 (lower_bound) if 0 is not the wanted result
        impl::cell_morton_idx(it->coord) == needle)
        return local_cell_index_type{
            std::distance(std::begin(_impl->m_p8est_cell_info), it)};
    else
        throw std::domain_error("Pos not in local domain.");
}

rank_type P4estGrid::position_to_rank(Vec3d pos)
{
    // Cell of pos might not be known on this process (not in
    // m_p8est_cell_info). Therefore, use the global first cell indices.
    auto it
        = std::upper_bound(std::begin(_impl->m_node_first_cell_idx),
                           std::end(_impl->m_node_first_cell_idx),
                           impl::pos_morton_idx(pos, _impl->m_inv_cell_size),
                           [](global_cell_index_type i,
                              global_cell_index_type idx) { return i < idx; });

    return std::distance(std::begin(_impl->m_node_first_cell_idx), it) - 1;
}

Vec3d P4estGrid::cell_size() const
{
    return _impl->m_cell_size;
}

Vec3i P4estGrid::grid_size() const
{
    return _impl->m_grid_size;
}

bool P4estGrid::repartition(CellMetric m,
                            CellCellMetric ccm,
                            Thunk exchange_start_callback)
{
    return _impl->repartition(m, ccm, exchange_start_callback);
}

global_cell_index_type
P4estGrid::global_hash(local_or_ghost_cell_index_type cellidx)
{
#ifdef GLOBAL_HASH_NEEDED
    return _impl->m_global_idx.at(_impl->locghost_to_p4est_index(cellidx));
#else
    return global_cell_index_type{0};
#endif
}

} // namespace grids
} // namespace repa
