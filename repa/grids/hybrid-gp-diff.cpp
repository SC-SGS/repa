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

//#ifdef HAVE_METIS

#include "hybrid-gp-diff.hpp"
#include "util/mpi_type.hpp"
#include <mpi.h>

namespace repa {
namespace grids {

HybridGPDiff::HybridGPDiff(const boost::mpi::communicator &comm,
                           Vec3d box_size,
                           double min_cell_size)
    : ParallelLCGrid(comm, box_size, min_cell_size),
      diff_impl(comm, box_size, min_cell_size),
      graph_impl(comm, box_size, min_cell_size),
      state(State::GRAPH),
      active_implementation(&graph_impl)
{
}

lidx HybridGPDiff::n_local_cells()
{
    return active_implementation->n_local_cells();
}

gidx HybridGPDiff::n_ghost_cells()
{
    return active_implementation->n_ghost_cells();
}

nidx HybridGPDiff::n_neighbors()
{
    return active_implementation->n_neighbors();
}

rank HybridGPDiff::neighbor_rank(nidx i)
{
    return active_implementation->neighbor_rank(i);
}

Vec3d HybridGPDiff::cell_size()
{
    return active_implementation->cell_size();
}

Vec3i HybridGPDiff::grid_size()
{
    return active_implementation->grid_size();
}

lgidx HybridGPDiff::cell_neighbor_index(lidx cellidx, int neigh)
{
    return active_implementation->cell_neighbor_index(cellidx, neigh);
}

std::vector<GhostExchangeDesc> HybridGPDiff::get_boundary_info()
{
    return active_implementation->get_boundary_info();
}

lidx HybridGPDiff::position_to_cell_index(double pos[3])
{
    return active_implementation->position_to_cell_index(pos);
}

rank HybridGPDiff::position_to_rank(double pos[3])
{
    return active_implementation->position_to_rank(pos);
}

nidx HybridGPDiff::position_to_neighidx(double pos[3])
{
    return active_implementation->position_to_neighidx(pos);
}

bool HybridGPDiff::repartition(CellMetric m,
                               CellCellMetric ccm,
                               Thunk exchange_start_callback)
{

    if (switch_to_state != state)
        switch_implementation();

    switch (state) {
    case State::GRAPH:
        if (comm_cart.rank() == 0)
            std::cout << "[Hybrid-GP-Diff] Repart GRAPH" << std::endl;
        break;
    case State::DIFF:
        if (comm_cart.rank() == 0)
            std::cout << "[Hybrid-GP-Diff] Repart DIFF" << std::endl;
        break;
    }

    return active_implementation->repartition(m, ccm, exchange_start_callback);
}

void HybridGPDiff::switch_implementation()
{
    if (state == switch_to_state)
        return;

    switch (switch_to_state) {
    case State::GRAPH:
        state = State::GRAPH;
        active_implementation = &graph_impl;
        std::copy(std::begin(diff_impl.partition),
                  std::end(diff_impl.partition),
                  std::begin(graph_impl.partition));
        MPI_Allreduce(MPI_IN_PLACE, graph_impl.partition.data(),
                      graph_impl.partition.size(),
                      MPI_ELEM_DECLTYPE_T(graph_impl.partition), MPI_MAX,
                      comm_cart);
        graph_impl.init();
        break;
    case State::DIFF:
        state = State::DIFF;
        active_implementation = &diff_impl;
        std::copy(std::begin(graph_impl.partition),
                  std::end(graph_impl.partition),
                  std::begin(diff_impl.partition));
        diff_impl.reinit();
        break;
    }
}

void HybridGPDiff::command(std::string s)
{
    // Only delegate the switch to later because if the user
    // were to switch the states not *directly* before a repartition
    // command, the grid would be inconsistent.
    // Because the active_implementation is changed to the implementation
    // that does not hold the current partitioning but something else.
    if (s == "graph" || s == "set graph") {
        switch_to_state = State::GRAPH;
    }
    else if (s == "diff" || s == "diffusion" || s == "set diff"
             || s == "set diffusion") {
        switch_to_state = State::DIFF;
    }
    else if (s == "toggle" || s == "switch") {
        switch_to_state = state == State::DIFF ? State::GRAPH : State::DIFF;
    }
}

} // namespace grids
} // namespace repa

//#endif // HAVE_METIS
