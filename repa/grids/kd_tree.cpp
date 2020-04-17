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

#include "kd_tree.hpp"
#include "util/linearize.hpp"
#include "util/neighbor_offsets.hpp"
#include "util/vec_arith.hpp"

#include <algorithm>
#include <cassert>
#include <kdpart/kdpart.h>
#include <numeric>

static int volume(const repa::Vec3i &v)
{
    return v[0] * v[1] * v[2];
}

static int volume(const repa::grids::Domain &dom)
{
    using namespace repa::util::vector_arithmetic;
    return volume(dom.second - dom.first);
}

static repa::Vec3i get_grid_dimensions(const repa::Vec3d &box_l,
                                       double max_range)
{
    using namespace repa::util::vector_arithmetic;
    if (max_range > ROUND_ERROR_PREC * box_l[0])
        return static_cast_vec<repa::Vec3i>(box_l / max_range);
    else
        return constant_vec3(1);
}

static repa::Vec3d get_cell_size(const repa::Vec3d &box_l,
                                 const repa::Vec3i &grid_dimensions)
{
    using namespace repa::util::vector_arithmetic;
    return box_l / grid_dimensions;
}

static bool does_domain_contain_cell(const repa::grids::Domain &domain,
                                     const repa::Vec3i &cell)
{
    using namespace repa::util::vector_arithmetic;
    return all(cell >= domain.first) && all(cell < domain.second);
}

static bool is_ghost_cell(const repa::Vec3i &cell,
                          const repa::Vec3i &ghostdomain_size)
{
    using namespace repa::util::vector_arithmetic;
    return any(cell == 0) || any(cell.as_expr() == ghostdomain_size - 1);
}

static repa::grids::Domain
to_ghost_domain_bounds(const repa::grids::Domain &domain)
{
    using namespace repa::util::vector_arithmetic;
    return {domain.first - 1, domain.second + 1};
}

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
static std::vector<repa::grids::Domain>
intersection_domains(const repa::grids::Domain &localdomain,
                     const repa::grids::Domain &ghostdomain,
                     const repa::Vec3i &global_grid_size,
                     bool ghostdomain_coords = false,
                     bool periodic_intersections_only = false)
{
    /*
#ifndef NDEBUG
    {
        // Preconditions:
        for (int dim = 0; dim < 3; dim++) {
            // Ensure valid size and position of localdomain
            assert(localdomain.first[dim] >= m_global_domain.first[dim]
                   && localdomain.first[dim] <= m_global_domain.second[dim]);
            assert(localdomain.second[dim] >= m_global_domain.first[dim]
                   && localdomain.second[dim] <= m_global_domain.second[dim]);
            assert(localdomain.first[dim] <= localdomain.second[dim]);
            // Ensure valid size and position of ghostdomain
            assert(ghostdomain.first[dim] >= m_global_ghostdomain.first[dim]
                   && ghostdomain.first[dim]
                          <= m_global_ghostdomain.second[dim]);
            assert(ghostdomain.second[dim] >= m_global_ghostdomain.first[dim]
                   && ghostdomain.second[dim]
                          <= m_global_ghostdomain.second[dim]);
            assert(ghostdomain.first[dim] <= ghostdomain.second[dim]);
        }
    }
#endif
    */
    using namespace repa::util::vector_arithmetic;

    // In perodic domains the local-domain must be shifted to perform
    // intersection tests with the ghostdomain across periodic domain bounds.
    // This is because the ghost-domain bounds can exceed the bounds of
    // the global domain.
    // In non-periodic domains only one check is necessary (=3^0).
    // In domains with one periodic dimension 3 checks are necessary (=3^1)
    // up to a maximum of 27 (=3^3) checks if alls dimensions are periodic.
    const repa::Vec3i periodicity{PERIODIC(0), PERIODIC(1), PERIODIC(2)};
    repa::grids::Domain neighborhood_to_check = {
        repa::Vec3i{0, 0, 0} - periodicity, repa::Vec3i{1, 1, 1} + periodicity};

    // Datastructure for gathering the results
    std::vector<repa::grids::Domain> intersection_domains;

    for (repa::rank_index_type nidx = 0; nidx < volume(neighborhood_to_check);
         nidx++) {
        // Determine neighbor offset
        repa::Vec3i neighbor_offset
            = repa::util::unlinearize(nidx, neighborhood_to_check.second
                                                - neighborhood_to_check.first)
              + neighborhood_to_check.first;

        // Check if the "default" intersection that isn't the result of shifting
        // the localdomain across periodic domain bounds should be excluded from
        // the result
        if (periodic_intersections_only && all(neighbor_offset == 0))
            continue;

        const repa::Vec3i num_cells_to_shift
            = neighbor_offset * global_grid_size;
        repa::grids::Domain shifted_localdomain
            = {localdomain.first + num_cells_to_shift,
               localdomain.second + num_cells_to_shift};

        if (any(shifted_localdomain.second <= ghostdomain.first)
            || any(ghostdomain.second <= shifted_localdomain.first))
            continue;

        // Calculate the actual intersection
        repa::grids::Domain intersection_domain;
        for (size_t dim = 0; dim < 3; ++dim) {
            intersection_domain.first[dim] = std::max(
                shifted_localdomain.first[dim], ghostdomain.first[dim]);
            intersection_domain.second[dim] = std::min(
                shifted_localdomain.second[dim], ghostdomain.second[dim]);
        }

        // If coords should be relative to the localdomain, they must be
        // shifted back
        if (!ghostdomain_coords) {
            intersection_domain.first -= num_cells_to_shift;
            intersection_domain.second -= num_cells_to_shift;
        }
        intersection_domains.push_back(intersection_domain);
    }
    return intersection_domains;
}

/**
 * Returns true if the given localdomain and the given ghostdomain
 * intersect. This includes intersections that are the result of periodic
 * domain bounds.
 */
static bool are_domains_intersecting(const repa::grids::Domain &localdomain,
                                     const repa::grids::Domain &ghostdomain,
                                     const repa::Vec3i &global_grid_size)
{
    return !intersection_domains(localdomain, ghostdomain, global_grid_size)
                .empty();
}

namespace repa {
namespace grids {

struct KDTreePrivateImpl {
    kdpart::PartTreeStorage t;

    kdpart::PartTreeStorage *operator->()
    {
        return &t;
    }

    KDTreePrivateImpl() = delete;
    KDTreePrivateImpl(const kdpart::PartTreeStorage &) = delete;
    KDTreePrivateImpl(kdpart::PartTreeStorage &&t) : t(t)
    {
    }
};

void KDTreeGrid::init_local_domain_bounds()
{
    using namespace util::vector_arithmetic;
    m_local_subdomain = (*m_kdtree)->subdomain_bounds(comm_cart.rank());
    m_local_ghostdomain = to_ghost_domain_bounds(m_local_subdomain);
    m_local_subdomain_size = m_local_subdomain.second - m_local_subdomain.first;
    m_local_ghostdomain_size
        = m_local_ghostdomain.second - m_local_ghostdomain.first;
}

void KDTreeGrid::init_nb_of_cells()
{
    m_nb_of_local_cells = volume(m_local_subdomain);
    m_nb_of_ghost_cells
        = volume(m_local_ghostdomain) - volume(m_local_subdomain);
}

void KDTreeGrid::init_index_permutations()
{
    // Number of local and ghost cells
    local_or_ghost_cell_index_type nb_of_total_cells
        = volume(m_local_ghostdomain);

    m_index_permutations.resize(nb_of_total_cells);
    m_index_permutations_inverse.resize(nb_of_total_cells);

    local_cell_index_type localidx = 0;
    local_or_ghost_cell_index_type ghostidx = m_nb_of_local_cells;
    for (local_or_ghost_cell_index_type cellidx = 0;
         cellidx < nb_of_total_cells; cellidx++) {
        const Vec3i cell = util::unlinearize(cellidx, m_local_ghostdomain_size);
        local_or_ghost_cell_index_type &idx
            = is_ghost_cell(cell, m_local_ghostdomain_size) ? ghostidx
                                                            : localidx;
        m_index_permutations_inverse[idx] = cellidx;
        m_index_permutations[cellidx] = idx++;
    }
}

// TODO return an iterator and not an elephant vector
std::vector<Vec3i> KDTreeGrid::cells(const std::vector<Domain> &domains)
{
    using namespace util::vector_arithmetic;
    std::vector<Vec3i> result;
    for (const Domain &domain : domains) {
        for (local_or_ghost_cell_index_type cellidx = 0;
             cellidx < volume(domain); cellidx++) {
            const Vec3i cell_coords
                = util::unlinearize(cellidx, domain.second - domain.first)
                  + domain.first;
            result.push_back(cell_coords);
        }
    }
    return result;
}

void KDTreeGrid::init_neighborhood_information()
{
    m_neighbor_processes_inverse.resize(comm_cart.size());
    std::fill(std::begin(m_neighbor_processes_inverse),
              std::end(m_neighbor_processes_inverse), UNKNOWN_RANK);

    (*m_kdtree)->walkp(
        [this](auto node) {
            return are_domains_intersecting({node.lu(), node.ro()},
                                            m_local_ghostdomain,
                                            m_global_domain_size);
        },
        // Explicitly using "this->" here for compatibility
        // with gcc 5
        [this](auto node) {
            if (!node.inner())
                this->init_neighborhood_information(node.rank());
        });
}

void KDTreeGrid::init_neighborhood_information(rank_type neighbor_rank)
{
    GhostExchangeDesc gexd;
    gexd.dest = neighbor_rank;

    // Initialize neighborhood information about ghostcells that the local
    // process is receiving from the neighbor process.
    const Domain neighbor_subdomain
        = (*m_kdtree)->subdomain_bounds(neighbor_rank);
    init_recv_cells(gexd, neighbor_subdomain);

    // Initialize neighborhood information about localcells that the local
    // neighbor process is sending to the neighbor process.
    const Domain neighbor_ghostdomain
        = to_ghost_domain_bounds(neighbor_subdomain);
    init_send_cells(gexd, neighbor_ghostdomain);

    if (!gexd.send.empty()) { // send.empty() == recv.empty(), so test only one.
        m_neighbor_processes_inverse[neighbor_rank]
            = m_neighbor_processes.size();
        m_neighbor_processes.push_back(neighbor_rank);
        m_boundary_info.emplace_back(std::move(gexd));
    }
    else {
        assert(neighbor_rank == comm_cart.rank());
    }
}

void KDTreeGrid::init_recv_cells(GhostExchangeDesc &gexd,
                                 const Domain &neighbor_subdomain)
{
    using namespace util::vector_arithmetic;
    // Get overlapping cells between the neighbor subdomain and the local
    // ghostdomain.
    std::vector<Vec3i> intersecting_cellvectors = cells(intersection_domains(
        neighbor_subdomain, m_local_ghostdomain, m_global_domain_size, true,
        gexd.dest == comm_cart.rank()));
    gexd.recv.reserve(intersecting_cellvectors.size());

    for (const Vec3i &intersecting_cellvector : intersecting_cellvectors) {
        const local_or_ghost_cell_index_type ghostidx
            = m_index_permutations[util::linearize(
                intersecting_cellvector - m_local_ghostdomain.first,
                m_local_ghostdomain_size)];

        // Convert cell-id to ghostcell-id
        assert(ghostidx >= m_nb_of_local_cells
               && ghostidx < m_nb_of_ghost_cells + m_nb_of_local_cells);

        // Update datastructures
        gexd.recv.push_back(ghostidx);
    }
}

void KDTreeGrid::init_send_cells(GhostExchangeDesc &gexd,
                                 const Domain &neighbor_ghostdomain)
{
    using namespace util::vector_arithmetic;
    // Get overlapping cells between the neighbor ghostdomain and the local
    // domain.
    std::vector<Vec3i> intersecting_cellvectors = cells(intersection_domains(
        m_local_subdomain, neighbor_ghostdomain, m_global_domain_size, false,
        gexd.dest == comm_cart.rank()));
    gexd.send.reserve(intersecting_cellvectors.size());

    for (const Vec3i &intersecting_cellvector : intersecting_cellvectors) {
        const local_cell_index_type cellidx
            = util::linearize(intersecting_cellvector - m_local_subdomain.first,
                              m_local_subdomain_size);

        // Update datastructure
        gexd.send.push_back(cellidx);
    }
}

void KDTreeGrid::clear_lookup_datastructures()
{
    m_neighbor_processes.clear();
    m_neighbor_processes_inverse.clear();
    m_index_permutations.clear();
    m_index_permutations_inverse.clear();
    m_boundary_info.clear();
}

void KDTreeGrid::reinitialize()
{
    clear_lookup_datastructures();
    init_local_domain_bounds();
    init_nb_of_cells();
    init_index_permutations();
    init_neighborhood_information();
}

KDTreeGrid::KDTreeGrid(const boost::mpi::communicator &comm,
                       Vec3d box_size,
                       double min_cell_size)
    : ParallelLCGrid(comm, box_size, min_cell_size),
      m_global_domain_size(get_grid_dimensions(box_l, max_range)),
      m_global_domain({Vec3i{0, 0, 0}, m_global_domain_size}),
      m_global_ghostdomain(to_ghost_domain_bounds(m_global_domain)),
      m_cell_size(get_cell_size(box_l, m_global_domain_size))
{
    // Use constant load to make initial tree evenly distributed
    auto load_function = [](Vec3i) { return 1; };
    m_kdtree = std::unique_ptr<KDTreePrivateImpl>(new KDTreePrivateImpl(
        kdpart::make_parttree(comm.size(), {0, 0, 0},
                              m_global_domain_size.as_array(), load_function,
                              kdpart::quality_splitting)));
    reinitialize();
}

local_cell_index_type KDTreeGrid::n_local_cells()
{
    return m_nb_of_local_cells;
}

ghost_cell_index_type KDTreeGrid::n_ghost_cells()
{
    return m_nb_of_ghost_cells;
}

rank_index_type KDTreeGrid::n_neighbors()
{
    return m_neighbor_processes.size();
}

rank_type KDTreeGrid::neighbor_rank(rank_index_type i)
{
    assert(i >= 0 && i < m_neighbor_processes.size());
    return m_neighbor_processes[i];
}

Vec3d KDTreeGrid::cell_size()
{
    return m_cell_size;
}

Vec3i KDTreeGrid::grid_size()
{
    return m_global_domain_size;
}

local_or_ghost_cell_index_type
KDTreeGrid::cell_neighbor_index(local_cell_index_type cellidx,
                                fs_neighidx neigh)
{
    // Precondition
    assert(cellidx >= 0 && cellidx < m_nb_of_local_cells);

    using namespace util::vector_arithmetic;
    // Get cellvector in local ghostdomain
    Vec3i neighbor_coord
        = util::unlinearize(m_index_permutations_inverse[cellidx],
                            m_local_ghostdomain_size)
          + util::NeighborOffsets3D::raw[neigh];

    // Linearization of vectors in ghostdomain require index permutations
    // to retain ordering requirements of lgidx indices.
    return m_index_permutations[util::linearize(neighbor_coord,
                                                m_local_ghostdomain_size)];
}

std::vector<GhostExchangeDesc> KDTreeGrid::get_boundary_info()
{
    return m_boundary_info;
}

local_cell_index_type KDTreeGrid::position_to_cell_index(Vec3d pos)
{
    using namespace util::vector_arithmetic;
    const Vec3i cell_coords = static_cast_vec<Vec3i>(pos / m_cell_size);

    if (!does_domain_contain_cell(m_local_subdomain, cell_coords)) {
        throw std::domain_error("Position not in local subdomain");
    }

    const local_cell_index_type cellidx = util::linearize(
        cell_coords - m_local_subdomain.first, m_local_subdomain_size);

    assert(cellidx >= 0 && cellidx < m_nb_of_local_cells);
    return cellidx;
}

rank_type KDTreeGrid::position_to_rank(Vec3d pos)
{
    using namespace util::vector_arithmetic;
    const Vec3i cell_coords = static_cast_vec<Vec3i>(pos / m_cell_size);
    return (*m_kdtree)->responsible_process(cell_coords.as_array());
}

rank_index_type KDTreeGrid::position_to_neighidx(Vec3d pos)
{
    using namespace util::vector_arithmetic;
    const Vec3i cell_coords = static_cast_vec<Vec3i>(pos / m_cell_size);
    const rank_type prank
        = (*m_kdtree)->responsible_process(cell_coords.as_array());
    const rank_index_type rank_idx = m_neighbor_processes_inverse[prank];

    if (rank_idx == UNKNOWN_RANK) {
        throw std::domain_error("Position not within neighbor a process");
    }

    return rank_idx;
}

bool KDTreeGrid::repartition(CellMetric m, CellCellMetric ccm, Thunk cb)
{
    UNUSED(ccm);
    const auto weights = m();
    assert(weights.size() == n_local_cells());

    m_kdtree = std::unique_ptr<KDTreePrivateImpl>(new KDTreePrivateImpl(
        kdpart::repart_parttree_par(m_kdtree->t, comm_cart, m())));
    cb();
    reinitialize();
    return true;
}

global_cell_index_type
KDTreeGrid::global_hash(local_or_ghost_cell_index_type cellidx)
{
    using namespace util::vector_arithmetic;
    // No need to define this away. Does currently not require extra data.
    const Vec3i idx3d = util::unlinearize(m_index_permutations_inverse[cellidx],
                                          m_local_ghostdomain_size);
    const Vec3i offset = m_local_ghostdomain.first - 1;
    const Vec3i gloidx3d = (idx3d + offset) % m_global_domain_size;
    return util::linearize(gloidx3d, m_global_domain_size);
}

} // namespace grids
} // namespace repa
