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
#include <boost/mpi.hpp>

#include "util/initial_partitioning.hpp"
#include "util/push_back_unique.hpp"

#ifndef NDEBUG
#define GLOMETHOD_DEBUG
#endif

namespace repa {
namespace grids {

local_cell_index_type GloMethod::n_local_cells() const
{
    return localCells;
}

ghost_cell_index_type GloMethod::n_ghost_cells() const
{
    return ghostCells;
}

rank_index_type GloMethod::n_neighbors() const
{
    return neighbors.size();
}

rank_type GloMethod::neighbor_rank(rank_index_type i) const
{
    return neighbors[i];
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
    assert(cellidx >= 0 && cellidx < n_local_cells());
    return global_to_local[gbox.neighbor(cells[cellidx], neigh)];
}

std::vector<GhostExchangeDesc> GloMethod::get_boundary_info()
{
    return exchangeVector;
}

local_cell_index_type GloMethod::position_to_cell_index(Vec3d pos)
{
    try {
        const auto c = global_to_local.at(gbox.cell_at_pos(pos));
        if (c >= n_local_cells())
            throw std::domain_error("Particle not in local subdomain");
        return c;
    }
    catch (const std::out_of_range &e) {
        throw std::domain_error("Particle not in local subdomain");
    }
}

rank_type GloMethod::position_to_rank(Vec3d pos)
{
    auto r = rank_of_cell(gbox.cell_at_pos(pos));

    if (r == UNKNOWN_RANK)
        throw std::runtime_error("Cell not in scope of process");
    else
        return r;
}

rank_index_type GloMethod::position_to_neighidx(Vec3d pos)
{
    rank_type rank = position_to_rank(pos);

    // Need to iterate neighbor_rank because GridBasedGrid provides a custom
    // implementation of it.
    for (rank_index_type i = 0; i < n_neighbors(); ++i) {
        if (neighbor_rank(i) == rank)
            return i;
    }
    throw std::domain_error("Position not within a neighbor process.");
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
    const global_cell_index_type nglocells = gbox.ncells();

    pre_init(firstcall);

    localCells = 0;
    ghostCells = 0;
    cells.clear();
    global_to_local.clear();
    neighbors.clear();

    // Extract the local cells from "partition".
    for (global_cell_index_type i = 0; i < nglocells; i++) {
        if (rank_of_cell(i) == comm_cart.rank()) {
            // Vector of own cells
            cells.push_back(i);
            // Index mapping from global to local
            global_to_local[i] = localCells;
            // Number of own cells
            localCells++;
        }
    }

    // Temporary storage for exchange descriptors.
    // Will be filled only for neighbors
    // and moved from later.
    std::vector<GhostExchangeDesc> tmp_ex_descs(comm_cart.size());

    // Determine ghost cells and communication volume
    for (local_cell_index_type i = 0; i < localCells; i++) {
        for (global_cell_index_type neighborIndex :
             gbox.full_shell_neigh_without_center(cells[i])) {
            const rank_type owner = rank_of_cell(neighborIndex);
            assert(owner != UNKNOWN_RANK);
            if (owner == comm_cart.rank())
                continue;

            init_new_foreign_cell(i, neighborIndex, owner);

            // Register "neighborIndex" as ghost cell
            if (global_to_local.find(neighborIndex)
                == std::end(global_to_local)) {
                cells.push_back(neighborIndex);
                global_to_local[neighborIndex] = localCells + ghostCells;
                ghostCells++;
            }
            else {
                // Must have been registered before as ghost cell to be received
                // from "owner".
                assert(std::find(std::begin(tmp_ex_descs[owner].recv),
                                 std::end(tmp_ex_descs[owner].recv),
                                 neighborIndex)
                       != std::end(tmp_ex_descs[owner].recv));
            }

            // Initialize exdesc and add "owner" as neighbor if unknown.
            if (tmp_ex_descs[owner].dest == -1) {
                neighbors.push_back(owner);
                tmp_ex_descs[owner].dest = owner;
            }

            util::push_back_unique(tmp_ex_descs[owner].recv, neighborIndex);
            util::push_back_unique(tmp_ex_descs[owner].send, cells[i]);
        }
    }

    // Move all existent exchange descriptors from "tmp_ex_descs" to
    // "exchangeVector".
    exchangeVector.clear();
    for (rank_type i = 0; i < comm_cart.size(); ++i) {
        if (tmp_ex_descs[i].dest != -1) {
            auto ed = std::move(tmp_ex_descs[i]);

            // Make sure, index ordering is the same on every process
            // and global to local index conversion
            std::sort(std::begin(ed.recv), std::end(ed.recv));
            std::transform(std::begin(ed.recv), std::end(ed.recv),
                           std::begin(ed.recv),
                           [this](global_cell_index_type i) {
                               return global_to_local[i];
                           });
            std::sort(std::begin(ed.send), std::end(ed.send));
            std::transform(std::begin(ed.send), std::end(ed.send),
                           std::begin(ed.send),
                           [this](global_cell_index_type i) {
                               return global_to_local[i];
                           });

            exchangeVector.push_back(std::move(ed));
        }
    }

#ifdef GLOMETHOD_DEBUG
    for (rank_type i = 0; i < comm_cart.size(); ++i) {
        if (tmp_ex_descs[i].dest != -1)
            assert(tmp_ex_descs[i].recv.size() == 0
                   && tmp_ex_descs[i].send.size() == 0);
    }
#endif

    post_init(firstcall);
}

global_cell_index_type
GloMethod::global_hash(local_or_ghost_cell_index_type cellidx)
{
    // No need to define this away. Does currently not require extra data.
    return cells[cellidx];
}

} // namespace grids
} // namespace repa
