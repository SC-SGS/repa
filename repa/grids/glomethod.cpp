/**
 * Copyright 2017-2019 Steffen Hirschmann
 * Copyright 2017-2018 Maximilian Wildbrett
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

#include "glomethod.hpp"
#include <algorithm>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/mpi.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>

#include "util/push_back_unique.hpp"
#include "util/range.hpp"

namespace repa {
namespace grids {

local_cell_index_type GloMethod::n_local_cells() const
{
    return local_cell_index_type{cell_store.local_cells().size()};
}

ghost_cell_index_type GloMethod::n_ghost_cells() const
{
    return ghost_cell_index_type{cell_store.ghost_cells().size()};
}

util::const_span<rank_type> GloMethod::neighbor_ranks() const
{
    return util::make_const_span(neighbors);
}

Vec3d GloMethod::cell_size() const
{
    return gbox.cell_size();
}

Vec3i GloMethod::grid_size() const
{
    return gbox.grid_size();
}

local_or_ghost_cell_index_type
GloMethod::cell_neighbor_index(local_cell_index_type cellidx, fs_neighidx neigh)
{
    return cell_store.as_local_or_ghost_index(
        gbox.neighbor(cell_store.as_global_index(cellidx), neigh));
}

util::const_span<GhostExchangeDesc> GloMethod::get_boundary_info()
{
    return util::make_const_span(exchangeVector);
}

local_cell_index_type GloMethod::position_to_cell_index(Vec3d pos)
{
    try {
        const auto c = cell_store.as_local_index(gbox.cell_at_pos(pos));
        if (!c)
            throw std::domain_error("Particle not in local subdomain");
        return *c;
    }
    catch (const std::out_of_range &e) {
        throw std::domain_error("Particle not in local subdomain");
    }
}

rank_type GloMethod::position_to_rank(Vec3d pos)
{
    return rank_of_cell(gbox.cell_at_pos(pos))
        .value_or_throw<std::runtime_error>("Cell not in scope of process");
}

/*
 * Repartition.
 * Every node is responsible for a certain range of cells along the
 * linearization defined by gbox. Every process evaluates the weights for its
 * cells and sends them to the process that is responsible for the graph node
 * that corresponds to the cell.
 * Partitioning is performed in parallel via ParMETIS.
 */
bool GloMethod::repartition(CellMetric m,
                            CellCellMetric ccm,
                            Thunk exchange_start_callback)
{
    bool result = sub_repartition(m, ccm);
    if (result) {
        // Position to rank is answered solely via the "partition" vector.
        // So the particle migration is ready to go.
        exchange_start_callback();
        init();
    }
    return result;
}

GloMethod::GloMethod(const boost::mpi::communicator &comm,
                     Vec3d box_size,
                     double min_cell_size,
                     ExtraParams ep)
    : ParallelLCGrid(comm, box_size, min_cell_size),
      gbox(box_size, min_cell_size),
      initial_partitioning(ep.init_part
                               ? util::parse_part_type(*ep.init_part)
                               : util::InitialPartitionType::CARTESIAN3D)
{
}

void GloMethod::after_construction()
{
    init(true);
}

GloMethod::~GloMethod()
{
}

/*
 * Rebuild the data structures describing subdomain and communication.
 */
void GloMethod::init(bool firstcall)
{
    pre_init(firstcall);
    cell_store.clear();
    neighbors.clear();

    // Extract the local cells from "partition".
    for (const auto i : gbox.global_cells()) {
        if (auto r = rank_of_cell(i); r && *r == comm_cart.rank()) {
            cell_store.push_back_local(i);
        }
    }

    //
    // In the following, we first collect the communication volume as global
    // cell indices. We do this to ensure the same sorting among all processes.
    // After sorting, we map these indices to local ones.
    //
    using GlobalCellIndices = std::vector<global_cell_index_type>;
    std::vector<GlobalCellIndices> recvvol_per_rank(comm_cart.size()),
        sendvol_per_rank(comm_cart.size());

    // Step 1. Determine ghost cells and communication volume
    for (const auto i : cell_store.local_cells()) {
        for (const global_cell_index_type neighborIndex :
             gbox.full_shell_neigh_without_center(
                 cell_store.as_global_index(i))) {
            const rank_type owner = rank_of_cell(neighborIndex).value();

            if (owner == comm_cart.rank())
                continue;

            init_new_foreign_cell(i, neighborIndex, owner);

            // Register "neighborIndex" as ghost cell
            if (!cell_store.holds_global_index(neighborIndex)) {
                cell_store.push_back_ghost(neighborIndex);
            }
            else {
                // Cannot be local cell, skipped above if this rank is owner.
                assert(cell_store.as_local_or_ghost_index(neighborIndex)
                           .is<ghost_cell_index_type>());
                // Must have been registered before as ghost cell to be received
                // from "owner".
                assert(std::find(std::begin(recvvol_per_rank[owner]),
                                 std::end(recvvol_per_rank[owner]),
                                 neighborIndex)
                       != std::end(recvvol_per_rank[owner]));
            }

            util::push_back_unique(neighbors, owner);

            util::push_back_unique(recvvol_per_rank[owner], neighborIndex);
            util::push_back_unique(sendvol_per_rank[owner],
                                   cell_store.as_global_index(i));
        }
    }

    // Step 2. Sort the indices globally uniquely and map them to local ones.
    exchangeVector.clear();
    for (rank_type i = 0; i < comm_cart.size(); ++i) {
        // Move it so that it is directly deleted at the end of this scope
        auto recv = std::move(recvvol_per_rank[i]);
        auto send = std::move(sendvol_per_rank[i]);

        assert(recv.empty() == send.empty());

        if (recv.empty())
            continue;

        exchangeVector.emplace_back(
            i,
            boost::copy_range<std::vector<ghost_cell_index_type>>(
                boost::sort(recv)
                | boost::adaptors::transformed(
                    cell_store.global_to_ghost_transformer())),
            boost::copy_range<std::vector<local_cell_index_type>>(
                boost::sort(send)
                | boost::adaptors::transformed(
                    cell_store.global_to_local_transformer())));
    }

    post_init(firstcall);
}

global_cell_index_type
GloMethod::global_hash(local_or_ghost_cell_index_type cellidx)
{
    // No need to define this away. Does currently not require extra data.
    return cell_store.as_global_index(cellidx);
}

} // namespace grids
} // namespace repa
