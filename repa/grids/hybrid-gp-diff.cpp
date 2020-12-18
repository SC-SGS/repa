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

#include "hybrid-gp-diff.hpp"
#include <boost/mpi/datatype.hpp>
#include <mpi.h>

#define MPI_ELEM_DECLTYPE_T(arr)                                               \
    boost::mpi::get_mpi_datatype(decltype(arr)::value_type{})

namespace repa {
namespace grids {

HybridGPDiff::HybridGPDiff(const boost::mpi::communicator &comm,
                           Vec3d box_size,
                           double min_cell_size,
                           ExtraParams ep)
    : ParallelLCGrid(comm, box_size, min_cell_size),
      diff_impl(comm, box_size, min_cell_size, ep),
      graph_impl(comm, box_size, min_cell_size, ep),
      state(State::GRAPH),
      switch_to_state(State::GRAPH),
      active_implementation(&graph_impl)
{
}

void HybridGPDiff::after_construction()
{
    active_implementation->after_construction();
}

local_cell_index_type HybridGPDiff::n_local_cells() const
{
    return active_implementation->n_local_cells();
}

ghost_cell_index_type HybridGPDiff::n_ghost_cells() const
{
    return active_implementation->n_ghost_cells();
}

util::const_span<rank_type> HybridGPDiff::neighbor_ranks() const
{
    return active_implementation->neighbor_ranks();
}

Vec3d HybridGPDiff::cell_size() const
{
    return active_implementation->cell_size();
}

Vec3i HybridGPDiff::grid_size() const
{
    return active_implementation->grid_size();
}

local_or_ghost_cell_index_type
HybridGPDiff::cell_neighbor_index(local_cell_index_type cellidx,
                                  fs_neighidx neigh)
{
    return active_implementation->cell_neighbor_index(cellidx, neigh);
}

util::const_span<GhostExchangeDesc> HybridGPDiff::get_boundary_info()
{
    return active_implementation->get_boundary_info();
}

local_cell_index_type HybridGPDiff::position_to_cell_index(Vec3d pos)
{
    return active_implementation->position_to_cell_index(pos);
}

rank_type HybridGPDiff::position_to_rank(Vec3d pos)
{
    return active_implementation->position_to_rank(pos);
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
        active_implementation = &graph_impl;
        // Copy mapping empty optionals to "-1".
        for (size_t i = 0; i < diff_impl.partition.size(); ++i)
            graph_impl.partition[i] = diff_impl.partition[i].value_or(-1);
        MPI_Allreduce(MPI_IN_PLACE, graph_impl.partition.data(),
                      graph_impl.partition.size(),
                      MPI_ELEM_DECLTYPE_T(graph_impl.partition), MPI_MAX,
                      comm_cart);
#ifndef NDEBUG
        // Check that all -1 are gone.
        for (const auto &el : graph_impl.partition)
            assert(el >= 0 && el < comm.size());
#endif
        graph_impl.init();
        break;
    case State::DIFF:
        active_implementation = &diff_impl;
        std::copy(std::begin(graph_impl.partition),
                  std::end(graph_impl.partition),
                  std::begin(diff_impl.partition));
#ifndef NDEBUG
        for (const auto &el : graph_impl.partition)
            assert(el >= 0 && el < comm.size());
#endif
        diff_impl.init();
        break;
    }

    state = switch_to_state;
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

    // Delegate commands prefixed with "diff:" or "graph:" to implementations
    if (s.substr(0, 5) == "diff:") {
        diff_impl.command(s.substr(5));
    }
    else if (s.substr(0, 6) == "graph:") {
        graph_impl.command(s.substr(6));
    }
}

global_cell_index_type
HybridGPDiff::global_hash(local_or_ghost_cell_index_type cellidx)
{
    return active_implementation->global_hash(cellidx);
}

} // namespace grids
} // namespace repa
