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

#include "util/ensure.hpp"
#include "util/push_back_unique.hpp"

#ifndef NDEBUG
#define GLOBOX_DEBUG
#endif

namespace repa {
namespace grids {

lidx GloMethod::n_local_cells()
{
    return localCells;
}

gidx GloMethod::n_ghost_cells()
{
    return ghostCells;
}

nidx GloMethod::n_neighbors()
{
    return neighbors.size();
}

rank GloMethod::neighbor_rank(nidx i)
{
    return neighbors[i];
}

Vec3d GloMethod::cell_size()
{
    return gbox.cell_size();
}

Vec3i GloMethod::grid_size()
{
    return gbox.grid_size();
}

lgidx GloMethod::cell_neighbor_index(lidx cellidx, fs_neighidx neigh)
{
    return global_to_local[gbox.neighbor(cells[cellidx], neigh)];
}

std::vector<GhostExchangeDesc> GloMethod::get_boundary_info()
{
    return exchangeVector;
}

lidx GloMethod::position_to_cell_index(Vec3d pos)
{
    if (position_to_rank(pos) != comm_cart.rank())
        throw std::domain_error("Particle not in local box");

    return global_to_local[gbox.cell_at_pos(pos)];
}

rank GloMethod::position_to_rank(Vec3d pos)
{
    auto r = partition[gbox.cell_at_pos(pos)];

    if (r < 0)
        throw std::runtime_error("Cell not in scope of process");
    else
        return r;
}

nidx GloMethod::position_to_neighidx(Vec3d pos)
{
    rank rank = position_to_rank(pos);
    auto ni = std::find(std::begin(neighbors), std::end(neighbors), rank);

    if (ni != std::end(neighbors))
        return std::distance(std::begin(neighbors), ni);
    else
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
             double min_cell_size)
    : ParallelLCGrid(comm, box_size, min_cell_size),
      gbox(box_size, min_cell_size)
{
    int nglocells = gbox.ncells();
    int ncells_per_proc = static_cast<int>(
        std::ceil(static_cast<double>(nglocells) / comm_cart.size()));

    // Initial partitioning
    partition.resize(nglocells);

    // Line-wise init
    for (int i = 0; i < nglocells; ++i) {
        partition[i] = i / ncells_per_proc;
    }

    // TODO: choose init?

    //// Init to equally sized boxes on Cartesian grid
    // int dims[3] = {0, 0, 0};
    // MPI_Dims_create(comm_cart.size(), 3, dims);

    // auto cellgrid = gbox.grid_size();
    // Vec3i cells_per_proc = {{
    //    static_cast<int>(std::ceil(static_cast<double>(cellgrid[0]) /
    //    dims[0])), static_cast<int>(std::ceil(static_cast<double>(cellgrid[1])
    //    / dims[1])),
    //    static_cast<int>(std::ceil(static_cast<double>(cellgrid[2]) /
    //    dims[2])),
    //}};

    // for (int i = 0; i < nglocells; ++i) {
    //  auto cellidx = gbox.unlinearize(i);
    //  // Transform cellidx to 3d proc coord
    //  for (int i = 0; i < 3; ++i)
    //    cellidx[i] /= cells_per_proc[i];
    //  int rank;
    //  MPI_Cart_rank(comm_cart, cellidx.data(), &rank);
    //  partition[i] = rank;
    //}
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
    const int nglocells = partition.size();

    pre_init(firstcall);

    localCells = 0;
    ghostCells = 0;
    cells.clear();
    global_to_local.clear();
    neighbors.clear();

    // Extract the local cells from "partition".
    for (int i = 0; i < nglocells; i++) {
        if (partition[i] == comm_cart.rank()) {
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
    for (int i = 0; i < localCells; i++) {
        for (int neighborIndex :
             gbox.full_shell_neigh_without_center(cells[i])) {
            rank owner = partition[neighborIndex];
            if (owner == comm_cart.rank())
                continue;

            init_new_foreign_cell(i, neighborIndex, owner);

            // Find ghost cells. Add only once to "cells" vector.
            if (global_to_local.find(neighborIndex)
                == std::end(global_to_local)) {
                // Add ghost cell to cells vector
                cells.push_back(neighborIndex);
                // Index mapping from global to ghost
                global_to_local[neighborIndex] = localCells + ghostCells;
                // Number of ghost cells
                ghostCells++;
            }

            // Initialize exdesc and add "rank" as neighbor if unknown.
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
    for (int i = 0; i < comm_cart.size(); ++i) {
        if (tmp_ex_descs[i].dest != -1) {
            auto ed = std::move(tmp_ex_descs[i]);

            // Make sure, index ordering is the same on every process
            // and global to local index conversion
            std::sort(std::begin(ed.recv), std::end(ed.recv));
            std::transform(std::begin(ed.recv), std::end(ed.recv),
                           std::begin(ed.recv),
                           [this](int i) { return global_to_local[i]; });
            std::sort(std::begin(ed.send), std::end(ed.send));
            std::transform(std::begin(ed.send), std::end(ed.send),
                           std::begin(ed.send),
                           [this](int i) { return global_to_local[i]; });

            exchangeVector.push_back(std::move(ed));
        }
    }

#ifdef GLOMETHOD_DEBUG
    for (int i = 0; i < comm_cart.size(); ++i) {
        if (tmp_ex_descs[i].dest != -1)
            ENSURE(tmp_ex_descs[i].recv.size() == 0
                   && tmp_ex_descs[i].send.size() == 0);
    }
#endif

    post_init(firstcall);
}

int GloMethod::global_hash(lgidx cellidx)
{
    // No need to define this away. Does currently not require extra data.
    return cells[cellidx];
}

} // namespace grids
} // namespace repa
