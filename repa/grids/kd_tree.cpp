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

//#ifdef HAVE_KDPART

#include "kd_tree.hpp"
#include "util/linearize.hpp"
#include "util/neighbor_offsets.hpp"
#include "util/vadd.hpp"

#include <algorithm>
#include <cassert>
#include <numeric>

namespace repa {
namespace grids {

Vec3i KDTreeGrid::grid_dimensions()
{
    Vec3i grid_dimensions = {1, 1, 1};
    for (int dim = 0; dim < 3; dim++) {
        if (max_range > ROUND_ERROR_PREC * box_l[dim]) {
            grid_dimensions[dim]
                = std::max<Vec3i::value_type>(box_l[dim] / max_range, 1);
        }
    }
    return grid_dimensions;
}

Vec3d KDTreeGrid::cell_dimensions(const Vec3i &grid_dimensions)
{
    Vec3d cell_dimensions;
    for (int dim = 0; dim < 3; dim++) {
        cell_dimensions[dim] = box_l[dim] / grid_dimensions[dim];
    }
    return cell_dimensions;
}

int KDTreeGrid::volume(Vec3i domain_size)
{
    return domain_size[0] * domain_size[1] * domain_size[2];
}

int KDTreeGrid::volume(Domain domain_bounds)
{
    return volume(domain_size(domain_bounds));
}

Domain KDTreeGrid::ghostdomain_bounds(const Domain &domain)
{
    Domain ghostdomain;
    for (int dim = 0; dim < 3; dim++) {
        if (domain.second[dim] > domain.first[dim]) {
            ghostdomain.first[dim] = domain.first[dim] - 1;
            ghostdomain.second[dim] = domain.second[dim] + 1;
        }
    }
    return ghostdomain;
}

Vec3i KDTreeGrid::domain_size(const Domain &domain)
{
    Vec3i domain_size;
    for (int dim = 0; dim < 3; dim++) {
        domain_size[dim] = domain.second[dim] - domain.first[dim];
    }
    return domain_size;
}

bool KDTreeGrid::is_ghost_cell(const Vec3i &cell, const Vec3i &ghostdomain_size)
{
    for (auto dim = 0; dim < 3; dim++) {
        if (cell[dim] == 0 || cell[dim] == ghostdomain_size[dim] - 1) {
            return true;
        }
    }
    return false;
}

Vec3i KDTreeGrid::absolute_position_to_cell_position(
    const Vec3d &absolute_position)
{
#ifndef NDEBUG
    for (int dim = 0; dim < 3; dim++) {
        assert(absolute_position[dim] >= 0);
    }
#endif
    Vec3i cell_position;
    for (int dim = 0; dim < 3; dim++) {
        cell_position[dim] = absolute_position[dim] / m_cell_dimensions[dim];
    }
    return cell_position;
}

void KDTreeGrid::init_local_domain_bounds()
{
    m_local_subdomain = m_kdtree.subdomain_bounds(comm_cart.rank());
    m_local_ghostdomain = ghostdomain_bounds(m_local_subdomain);
    m_local_subdomain_size = domain_size(m_local_subdomain);
    m_local_ghostdomain_size = domain_size(m_local_ghostdomain);
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
    lgidx nb_of_total_cells = volume(m_local_ghostdomain);

    m_index_permutations.resize(nb_of_total_cells);
    m_index_permutations_inverse.resize(nb_of_total_cells);

    lidx localidx = 0;
    lgidx ghostidx = m_nb_of_local_cells;
    for (lgidx cellidx = 0; cellidx < nb_of_total_cells; cellidx++) {
        Vec3i cell = util::unlinearize(cellidx, m_local_ghostdomain_size);
        if (is_ghost_cell(cell, m_local_ghostdomain_size)) {
            m_index_permutations_inverse[ghostidx] = cellidx;
            m_index_permutations[cellidx] = ghostidx++;
        }
        else {
            m_index_permutations_inverse[localidx] = cellidx;
            m_index_permutations[cellidx] = localidx++;
        }
    }
}

std::vector<Domain>
KDTreeGrid::intersection_domains(const Domain &localdomain,
                                 const Domain &ghostdomain,
                                 bool ghostdomain_coords,
                                 bool periodic_intersections_only) const
{
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

    // In perodic domains the local-domain must be shifted to perform
    // intersection tests with the ghostdomain across periodic domain bounds.
    // This is because the ghost-domain bounds can exceed the bounds of
    // the global domain.
    // In non-periodic domains only one check is necessary (=3^0).
    // In domains with one periodic dimension 3 checks are necessary (=3^1)
    // up to a maximum of 27 (=3^3) checks if alls dimensions are periodic.
    Domain neighborhood_to_check = {{0, 0, 0}, {1, 1, 1}};
    for (int dim = 0; dim < 3; dim++) {
        if (PERIODIC(dim)) {
            neighborhood_to_check.first[dim]--;
            neighborhood_to_check.second[dim]++;
        }
    }

    // Datastructure for gathering the results
    std::vector<Domain> intersection_domains;

    for (rank_index_type nidx = 0; nidx < volume(neighborhood_to_check); nidx++) {
        // Determine neighbor offset
        Vec3i neighbor_offset
            = util::unlinearize(nidx, domain_size(neighborhood_to_check));
        for (int dim = 0; dim < 3; dim++) {
            neighbor_offset[dim] += neighborhood_to_check.first[dim];
        }

        // Check if the "default" intersection that isn't the result of shifting
        // the localdomain across periodic domain bounds should be excluded from
        // the result
        if (periodic_intersections_only && neighbor_offset == Vec3i{0, 0, 0}) {
            continue;
        }

        Domain intersection_domain;
        Domain shifted_localdomain;
        for (int dim = 0; dim < 3; dim++) {
            // Shift local domain along current dimension
            int num_cells_to_shift
                = neighbor_offset[dim] * m_global_domain_size[dim];
            shifted_localdomain.first[dim]
                = localdomain.first[dim] + num_cells_to_shift;
            shifted_localdomain.second[dim]
                = localdomain.second[dim] + num_cells_to_shift;

            // Check for intersection on current dimension
            if (shifted_localdomain.second[dim] <= ghostdomain.first[dim]
                || ghostdomain.second[dim] <= shifted_localdomain.first[dim]) {
                goto no_intersection; // ~ continue-statement on outer loop
            }

            // Calculate bounds of the intersection domain on the
            // current dimension.
            intersection_domain.first[dim] = std::max(
                shifted_localdomain.first[dim], ghostdomain.first[dim]);
            intersection_domain.second[dim] = std::min(
                shifted_localdomain.second[dim], ghostdomain.second[dim]);

            // If coords should be relative to the localdomain, they must be
            // shifted back
            if (!ghostdomain_coords) {
                intersection_domain.first[dim] -= num_cells_to_shift;
                intersection_domain.second[dim] -= num_cells_to_shift;
            }
        }

        intersection_domains.push_back(intersection_domain);
    no_intersection:;
    }

    return intersection_domains;
}

bool KDTreeGrid::are_domains_intersecting(const Domain &localdomain,
                                          const Domain &ghostdomain) const
{
    return intersection_domains(localdomain, ghostdomain).empty() == false;
}

bool KDTreeGrid::domain_contains_cell(const Domain &domain, const Vec3i &cell)
{
    for (int dim = 0; dim < 3; dim++) {
        if (cell[dim] < domain.first[dim] || cell[dim] >= domain.second[dim]) {
            return false;
        }
    }
    return true;
}

// TODO return an iterator and not an elephant vector
std::vector<Vec3i> KDTreeGrid::cells(const std::vector<Domain> &domains)
{
    std::vector<Vec3i> result;
    for (const Domain &domain : domains) {
        for (lgidx cellidx = 0; cellidx < volume(domain); cellidx++) {
            Vec3i cellvector = util::unlinearize(cellidx, domain_size(domain));
            for (int dim = 0; dim < 3; dim++) {
                cellvector[dim] += domain.first[dim];
            }
            result.push_back(cellvector);
        }
    }
    return result;
}

void KDTreeGrid::init_neighborhood_information()
{
    m_neighbor_processes_inverse.resize(comm_cart.size(), -1);

    m_kdtree.walkp(
        // Explicitly using "this->" here for compatibility
        // with gcc 5
        [this](auto node) {
            return this->are_domains_intersecting({node.lu(), node.ro()},
                                                  m_local_ghostdomain);
        },
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
    Domain neighbor_subdomain = m_kdtree.subdomain_bounds(neighbor_rank);
    init_recv_cells(gexd, neighbor_subdomain);

    // Initialize neighborhood information about localcells that the local
    // neighbor process is sending to the neighbor process.
    Domain neighbor_ghostdomain = ghostdomain_bounds(neighbor_subdomain);
    init_send_cells(gexd, neighbor_ghostdomain);

    if (gexd.send.empty()) { // send.empty() == recv.empty(), so test only one.
        // This can only be the case if neighbor_rank == comm_cart.rank()
        // and there are enough processes such that no self communication
        // is necessary.
        if (neighbor_rank != comm_cart.rank()) {
            throw std::runtime_error("Empty ghost exchange detected");
        }

        // Simply do not add this descriptor to the neighbors
        return;
    }
    // Add process to neighbor processes
    m_neighbor_processes_inverse[neighbor_rank] = m_neighbor_processes.size();
    m_neighbor_processes.push_back(neighbor_rank);
    m_boundary_info.emplace_back(std::move(gexd));
}

void KDTreeGrid::init_recv_cells(GhostExchangeDesc &gexd,
                                 const Domain &neighbor_subdomain)
{
    // Get overlapping cells between the neighbor subdomain and the local
    // ghostdomain.
    std::vector<Vec3i> intersecting_cellvectors
        = cells(intersection_domains(neighbor_subdomain, m_local_ghostdomain,
                                     true, gexd.dest == comm_cart.rank()));
    gexd.recv.reserve(intersecting_cellvectors.size());

    for (const Vec3i &intersecting_cellvector : intersecting_cellvectors) {
        // Convert global cellvector to local cellvector relative to
        // ghostdomain.
        Vec3i local_cellvector = intersecting_cellvector;
        for (int dim = 0; dim < 3; dim++) {
            local_cellvector[dim] -= m_local_ghostdomain.first[dim];
        }

        // Convert local cellvector to local/ghostcell-id
        lgidx ghostidx = m_index_permutations[util::linearize(
            local_cellvector, m_local_ghostdomain_size)];

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
    // Get overlapping cells between the neighbor ghostdomain and the local
    // domain.
    std::vector<Vec3i> intersecting_cellvectors
        = cells(intersection_domains(m_local_subdomain, neighbor_ghostdomain,
                                     false, gexd.dest == comm_cart.rank()));
    gexd.send.reserve(intersecting_cellvectors.size());

    for (const Vec3i &intersecting_cellvector : intersecting_cellvectors) {
        // Convert global cellvector to local cellvector relative to
        // ghostdomain.
        Vec3i local_cellvector = intersecting_cellvector;
        for (int dim = 0; dim < 3; dim++) {
            local_cellvector[dim] -= m_local_subdomain.first[dim];
        }

        // Convert local cellvector to local cell-id
        lidx cellidx
            = util::linearize(local_cellvector, m_local_subdomain_size);

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
      m_global_domain_size(grid_dimensions()),
      m_global_domain({{0, 0, 0}, m_global_domain_size}),
      m_global_ghostdomain(ghostdomain_bounds(m_global_domain)),
      m_global_ghostdomain_size(domain_size(m_global_ghostdomain)),
      m_cell_dimensions(cell_dimensions(m_global_domain_size))
{
    int nb_of_subdomains = comm_cart.size();

    // Use constant load to make initial tree evenly distributed
    auto load_function = [](Vec3i) { return 1; };

    m_kdtree = kdpart::make_parttree(nb_of_subdomains, {0, 0, 0},
                                     m_global_domain_size.as_array(),
                                     load_function, kdpart::quality_splitting);

    reinitialize();
}

lidx KDTreeGrid::n_local_cells()
{
    return m_nb_of_local_cells;
}

gidx KDTreeGrid::n_ghost_cells()
{
    return m_nb_of_ghost_cells;
}

rank_index_type KDTreeGrid::n_neighbors()
{
    return m_neighbor_processes.size();
}

rank_type KDTreeGrid::neighbor_rank(rank_index_type i)
{
    if (i < 0 || i >= m_neighbor_processes.size()) {
        throw std::domain_error("Invalid neighbor index.");
    }
    return m_neighbor_processes[i];
}

Vec3d KDTreeGrid::cell_size()
{
    return m_cell_dimensions;
}

Vec3i KDTreeGrid::grid_size()
{
    return m_global_domain_size;
}

lgidx KDTreeGrid::cell_neighbor_index(lidx cellidx, fs_neighidx neigh)
{
    // Preconditions
    if (cellidx < 0 || cellidx >= m_nb_of_local_cells)
        throw std::domain_error("Invalid cell index");

    // Get cellvector in local ghostdomain
    Vec3i cellvector = util::unlinearize(m_index_permutations_inverse[cellidx],
                                         m_local_ghostdomain_size);

    // Shift cellvector to the choosen neighbor
    const Vec3i &offset = util::NeighborOffsets3D::raw[neigh];
    for (int dim = 0; dim < 3; dim++) {
        cellvector[dim] += offset[dim];
    }

    // Linearization of vectors in ghostdomain require index permutations
    // to retain ordering requirements of lgidx indices.
    return m_index_permutations[util::linearize(cellvector,
                                                m_local_ghostdomain_size)];
}

std::vector<GhostExchangeDesc> KDTreeGrid::get_boundary_info()
{
    return m_boundary_info;
}

lidx KDTreeGrid::position_to_cell_index(Vec3d pos)
{
    Vec3i cellvector = absolute_position_to_cell_position(pos);

    // Precondition: Position must be within local subdomain
    if (!domain_contains_cell(m_local_subdomain, cellvector)) {
        throw std::domain_error("Position not in local subdomain");
    }

    // Make cellvector relative to subdomain
    for (int dim = 0; dim < 3; dim++) {
        cellvector[dim] -= m_local_subdomain.first[dim];
    }

    // Linearize cell vector
    lidx cellidx = util::linearize(cellvector, m_local_subdomain_size);

    return cellidx;
}

rank_type KDTreeGrid::position_to_rank(Vec3d pos)
{
    Vec3i cellvector = absolute_position_to_cell_position(pos);

    // Precondition: Position must be within global domain
    if (!domain_contains_cell(m_global_domain, cellvector)) {
        throw std::domain_error("Position not within global domain");
    }

    return m_kdtree.responsible_process(cellvector.as_array());
}

rank_index_type KDTreeGrid::position_to_neighidx(Vec3d pos)
{
    Vec3i cellvector = absolute_position_to_cell_position(pos);

    // Precondition
    if (!domain_contains_cell(m_global_ghostdomain, cellvector)) {
        throw std::domain_error("Position is not in the global ghostdomain");
    }

    rank_type prank = m_kdtree.responsible_process(cellvector.as_array());
    rank_index_type rank_idx = m_neighbor_processes_inverse[prank];
    if (rank_idx == -1) {
        throw std::domain_error("Position not within neighbor a process");
    }
    return rank_idx;
}

bool KDTreeGrid::repartition(CellMetric m, CellCellMetric ccm, Thunk cb)
{
    UNUSED(ccm);
    auto weights = m();
    if (weights.size() != n_local_cells()) {
        throw std::runtime_error(
            "Metric only supplied " + std::to_string(weights.size())
            + "weights. Necessary: " + std::to_string(n_local_cells()));
    }

    m_kdtree = repart_parttree_par(m_kdtree, comm_cart, m());
    cb();
    reinitialize();
    return true;
}

gloidx KDTreeGrid::global_hash(lgidx cellidx)
{
    // No need to define this away. Does currently not require extra data.
    Vec3i idx3d = util::unlinearize(m_index_permutations_inverse[cellidx],
                                    m_local_ghostdomain_size);
    Vec3i offset;
    for (size_t i = 0; i < idx3d.size(); ++i)
        offset[i] = m_local_subdomain.first[i] - 1; // -1 for ghost

    auto gloidx3d = util::vadd_mod(idx3d, offset, m_global_domain_size);
    return util::linearize(gloidx3d, m_global_domain_size);
}

} // namespace grids
} // namespace repa

//#endif // HAVE_KDPART
