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

#include <algorithm>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext.hpp>
#include <boost/range/numeric.hpp>
#include <cassert>
#include <numeric>

#include <kdpart/kdpart.h>

#include "kd_tree.hpp"
#include "util/linearize.hpp"
#include "util/neighbor_offsets.hpp"
#include "util/range.hpp"
#include "util/vec_arith.hpp"

namespace repa {
namespace grids {

namespace { // Anonymous namespace for internal linkage

int volume(const KDTreeGrid::Domain &dom)
{
    using namespace repa::util::vector_arithmetic;
    return product(dom.second - dom.first);
}

KDTreeGrid::Domain to_ghost_domain_bounds(const KDTreeGrid::Domain &domain)
{
    using namespace util::vector_arithmetic;
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
std::vector<KDTreeGrid::Domain>
intersection_domains(const KDTreeGrid::Domain &localdomain,
                     const KDTreeGrid::Domain &ghostdomain,
                     const Vec3i &global_grid_size,
                     bool ghostdomain_coords = false,
                     bool periodic_intersections_only = false)
{
    using namespace util::vector_arithmetic;

    // In perodic domains the local-domain must be shifted to perform
    // intersection tests with the ghostdomain across periodic domain bounds.
    // This is because the ghost-domain bounds can exceed the bounds of
    // the global domain.
    // In non-periodic domains only one check is necessary (=3^0).
    // In domains with one periodic dimension 3 checks are necessary (=3^1)
    // up to a maximum of 27 (=3^3) checks if alls dimensions are periodic.
    const Vec3i periodicity{PERIODIC(0), PERIODIC(1), PERIODIC(2)};
    KDTreeGrid::Domain neighborhood_to_check = {
        repa::Vec3i{0, 0, 0} - periodicity, repa::Vec3i{1, 1, 1} + periodicity};

    // Datastructure for gathering the results
    std::vector<KDTreeGrid::Domain> intersection_domains;

    for (const int nidx : boost::irange(volume(neighborhood_to_check))) {
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
        KDTreeGrid::Domain shifted_localdomain
            = {localdomain.first + num_cells_to_shift,
               localdomain.second + num_cells_to_shift};

        if (any(shifted_localdomain.second <= ghostdomain.first)
            || any(ghostdomain.second <= shifted_localdomain.first))
            continue;

        // Calculate the actual intersection
        KDTreeGrid::Domain intersection_domain;
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
bool are_domains_intersecting(const KDTreeGrid::Domain &localdomain,
                              const KDTreeGrid::Domain &ghostdomain,
                              const Vec3i &global_grid_size)
{
    return !intersection_domains(localdomain, ghostdomain, global_grid_size)
                .empty();
}

/** Appends all global 3d cell indices in the range of "domain" to vector
 * "result". */
std::vector<global_cell_index_type>
global_cellindices_of_domain(const KDTreeGrid::Domain &domain,
                             std::vector<global_cell_index_type> &result,
                             const util::box_global_index_storage &cell_store)
{
    using namespace util::vector_arithmetic;
    const auto ncells = volume(domain);
    for (int cell_no = 0; cell_no < ncells; ++cell_no) {
        const Vec3i cell_coords
            = util::unlinearize(cell_no, domain.second - domain.first)
              + domain.first;
        result.push_back(cell_store.get_canonical_representant(cell_coords));
    }
    return result;
}

/** Returns all global 3d cell indices in the range of the domains. */
std::vector<global_cell_index_type>
global_cellindices_of_domains(const std::vector<KDTreeGrid::Domain> &domains,
                              const util::box_global_index_storage &cell_store)
{
    std::vector<global_cell_index_type> result;
    // Expected size of all domains is their volume.
    const auto siz
        = boost::accumulate(domains | boost::adaptors::transformed(volume), 0);
    result.reserve(siz);
    for (const auto &domain : domains) {
        global_cellindices_of_domain(domain, result, cell_store);
    }

    // Sort and kill duplicates
    // Sorting is important: Ensures same ordering of cells on all processes.
    boost::erase(result,
                 boost::unique<boost::return_found_end>(boost::sort(result)));
    return result;
}

} // namespace

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

void KDTreeGrid::init_neighborhood_information()
{
    (*m_kdtree)->walkp(
        [this](auto node) {
            return are_domains_intersecting(
                {node.lu(), node.ro()}, m_local_ghostdomain, gbox.grid_size());
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
    // No self-communication
    if (neighbor_rank == comm.rank())
        return;

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
    // Get overlapping cells between the neighbor subdomain and the local
    // ghostdomain.
    const auto cellindices = global_cellindices_of_domains(
        intersection_domains(neighbor_subdomain, m_local_ghostdomain,
                             gbox.grid_size(), true,
                             gexd.dest == comm_cart.rank()),
        cell_store);

    gexd.recv.reserve(cellindices.size());
    // Leave out local cells
    for (const global_cell_index_type &gloidx : cellindices)
        if (const auto ghostidx = cell_store.as_ghost_index(gloidx))
            gexd.recv.push_back(*ghostidx);
}

void KDTreeGrid::init_send_cells(GhostExchangeDesc &gexd,
                                 const Domain &neighbor_ghostdomain)
{
    // Get overlapping cells between the neighbor ghostdomain and the local
    // domain.
    const auto cellindices = global_cellindices_of_domains(
        intersection_domains(m_local_subdomain, neighbor_ghostdomain,
                             gbox.grid_size(), false,
                             gexd.dest == comm_cart.rank()),
        cell_store);

    gexd.send.reserve(cellindices.size());
    for (const global_cell_index_type &gloidx : cellindices)
        gexd.send.push_back(cell_store.as_local_index(gloidx).ensure_value());
}

void KDTreeGrid::reinitialize()
{
    m_neighbor_processes.clear();
    m_boundary_info.clear();

    using namespace util::vector_arithmetic;
    m_local_subdomain = (*m_kdtree)->subdomain_bounds(comm_cart.rank());
    m_local_ghostdomain = to_ghost_domain_bounds(m_local_subdomain);

    cell_store.clear();
    cell_store.init(gbox.grid_size(),
                    m_local_subdomain.second - m_local_subdomain.first,
                    m_local_subdomain.first);

    init_neighborhood_information();
}

KDTreeGrid::KDTreeGrid(const boost::mpi::communicator &comm,
                       Vec3d box_size,
                       double min_cell_size)
    : ParallelLCGrid(comm, box_size, min_cell_size),
      gbox(box_size, min_cell_size)
{
    // Use constant load to make initial tree evenly distributed
    auto load_function = [](Vec3i) { return 1; };
    m_kdtree = std::make_unique<KDTreePrivateImpl>(kdpart::make_parttree(
        comm.size(), {0, 0, 0}, gbox.grid_size().as_array(), load_function,
        kdpart::quality_splitting));
    reinitialize();
}

local_cell_index_type KDTreeGrid::n_local_cells() const
{
    return local_cell_index_type{cell_store.local_cells().size()};
}

ghost_cell_index_type KDTreeGrid::n_ghost_cells() const
{
    return ghost_cell_index_type{cell_store.ghost_cells().size()};
}

util::const_span<rank_type> KDTreeGrid::neighbor_ranks() const
{
    return util::make_const_span(m_neighbor_processes);
}

Vec3d KDTreeGrid::cell_size() const
{
    return gbox.cell_size();
}

Vec3i KDTreeGrid::grid_size() const
{
    return gbox.grid_size();
}

local_or_ghost_cell_index_type
KDTreeGrid::cell_neighbor_index(local_cell_index_type cellidx,
                                fs_neighidx neigh)
{
    const auto gloidx = cell_store.as_global_index(cellidx);
    const auto neighidx = gbox.neighbor(gloidx, neigh);
    return cell_store.as_local_or_ghost_index(neighidx);
}

util::const_span<GhostExchangeDesc> KDTreeGrid::get_boundary_info()
{
    return util::make_const_span(m_boundary_info);
}

local_cell_index_type KDTreeGrid::position_to_cell_index(Vec3d pos)
{
    using namespace util::vector_arithmetic;
    const Vec3i cell_coords = static_cast_vec<Vec3i>(pos / gbox.cell_size());

    if (!(all(cell_coords >= m_local_subdomain.first)
          && all(cell_coords < m_local_subdomain.second)))
        throw std::domain_error("Position not in local subdomain");

    if (const auto idx = cell_store.as_local_index(cell_coords)) {
        assert(position_to_rank(pos) == comm.rank());
        return *idx;
    }
    else {
        // Not reached. Checked above that "pos" is in the local subdomain
        ensure_not_reached();
    }
}

rank_type KDTreeGrid::position_to_rank(Vec3d pos)
{
    using namespace util::vector_arithmetic;
    const Vec3i cell_coords = static_cast_vec<Vec3i>(pos / gbox.cell_size());
    return (*m_kdtree)->responsible_process(cell_coords.as_array());
}

bool KDTreeGrid::repartition(CellMetric m, CellCellMetric ccm, Thunk cb)
{
    UNUSED(ccm);
    const auto weights = m();
    assert(weights.size() == local_cells().size());

    m_kdtree = std::make_unique<KDTreePrivateImpl>(
        kdpart::repart_parttree_par(m_kdtree->t, comm_cart, m()));
    cb();
    reinitialize();
    return true;
}

global_cell_index_type
KDTreeGrid::global_hash(local_or_ghost_cell_index_type cellidx)
{
    return cell_store.as_global_index(cellidx);
}

} // namespace grids
} // namespace repa
