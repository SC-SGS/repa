
#ifdef HAVE_METIS

#include <mpi.h>
#include "hybrid-gp-diff.hpp"

namespace generic_dd {
namespace grids {

HybridGPDiff::HybridGPDiff()
    : state(State::GRAPH), active_implementation(&graph_impl) {}

lidx HybridGPDiff::n_local_cells() {
  return active_implementation->n_local_cells();
}

gidx HybridGPDiff::n_ghost_cells() {
  return active_implementation->n_ghost_cells();
}

nidx HybridGPDiff::n_neighbors() { return active_implementation->n_neighbors(); }

rank HybridGPDiff::neighbor_rank(nidx i) {
  return active_implementation->neighbor_rank(i);
}

std::array<double, 3> HybridGPDiff::cell_size() {
  return active_implementation->cell_size();
}

std::array<int, 3> HybridGPDiff::grid_size() {
  return active_implementation->grid_size();
}

lgidx HybridGPDiff::cell_neighbor_index(lidx cellidx, int neigh) {
  return active_implementation->cell_neighbor_index(cellidx, neigh);
}

std::vector<GhostExchangeDesc> HybridGPDiff::get_boundary_info() {
  return active_implementation->get_boundary_info();
}

lidx HybridGPDiff::position_to_cell_index(double pos[3]) {
  return active_implementation->position_to_cell_index(pos);
}

rank HybridGPDiff::position_to_rank(double pos[3]) {
  return active_implementation->position_to_rank(pos);
}

nidx HybridGPDiff::position_to_neighidx(double pos[3]) {
  return active_implementation->position_to_neighidx(pos);
}

bool HybridGPDiff::repartition(const repart::Metric &m,
                               std::function<void()> exchange_start_callback) {

  if (switch_to_state != state)
    switch_implementation();

  switch (state) {
  case State::GRAPH:
    if (this_node == 0)
      std::cout << "[Hybrid-GP-Diff] Repart GRAPH" << std::endl;
    break;
  case State::DIFF:
    if (this_node == 0)
      std::cout << "[Hybrid-GP-Diff] Repart DIFF" << std::endl;
    break;
  }

  return active_implementation->repartition(m, exchange_start_callback);
}


template <typename T>
struct mpi_type {};

template <>
struct mpi_type<int32_t> {
  static const MPI_Datatype type;
};
const MPI_Datatype mpi_type<int32_t>::type = MPI_INT32_T;

template <>
struct mpi_type<int64_t> {
  static const MPI_Datatype type;
};
const MPI_Datatype mpi_type<int64_t>::type = MPI_INT64_T;


void HybridGPDiff::switch_implementation() {
  if (state == switch_to_state)
    return;

  switch (switch_to_state) {
  case State::GRAPH:
    state = State::GRAPH;
    active_implementation = &graph_impl;
    std::copy(std::begin(diff_impl.partition), std::end(diff_impl.partition),
              std::begin(graph_impl.partition));
    MPI_Allreduce(MPI_IN_PLACE, graph_impl.partition.data(),
                  graph_impl.partition.size(), mpi_type<decltype(graph_impl.partition)::value_type>::type, MPI_MAX, comm_cart);
    graph_impl.init();
    break;
  case State::DIFF:
    state = State::DIFF;
    active_implementation = &diff_impl;
    std::copy(std::begin(graph_impl.partition), std::end(graph_impl.partition),
              std::begin(diff_impl.partition));
    diff_impl.reinit();
    break;
  }
}

void HybridGPDiff::command(std::string s) {
  // Only delegate the switch to later because if the user
  // were to switch the states not *directly* before a repartition
  // command, the grid would be inconsistent.
  // Because the active_implementation is changed to the implementation
  // that does not hold the current partitioning but something else.
  if (s == "graph" || s == "set graph") {
    switch_to_state = State::GRAPH;
  } else if (s == "diff" || s == "diffusion" || s == "set diff" ||
             s == "set diffusion") {
    switch_to_state = State::DIFF;
  } else if (s == "toggle" || s == "switch") {
    switch_to_state = state == State::DIFF ? State::GRAPH : State::DIFF;
  }
}

} // namespace grids
} // namespace generic_dd

#endif // HAVE_METIS
