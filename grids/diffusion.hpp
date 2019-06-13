#ifndef _GRIDS_DIFF_HPP
#define _GRIDS_DIFF_HPP

#include "../pargrid.hpp"
#include "globox.hpp"
#include <array>
#include <unordered_map>
#include <mpi.h>
#include <vector>

namespace generic_dd {
namespace grids {

struct Diffusion : public ParallelLCGrid {

  Diffusion();
  lidx n_local_cells() override;
  gidx n_ghost_cells() override;
  nidx n_neighbors() override;
  rank neighbor_rank(nidx i) override;
  std::array<double, 3> cell_size() override;
  std::array<int, 3> grid_size() override;
  lgidx cell_neighbor_index(lidx cellidx, int neigh) override;
  std::vector<GhostExchangeDesc> get_boundary_info() override;
  lidx position_to_cell_index(double pos[3]) override;
  rank position_to_rank(double pos[3]) override;
  nidx position_to_neighidx(double pos[3]) override;
  bool repartition(const repart::Metric &m,
                   std::function<void()> exchange_start_callback) override;

  struct NeighSend {
    int basecell;
    std::array<int, 26> neighranks;
  };
private:
  friend struct HybridGPDiff; // Needs access to "partition" vector
  // Number of local cells of this process(partition)
  lidx localCells;
  // Number of ghost cells of this process(partition)
  int ghostCells;
  // Vector with all neighbour ranks
  std::vector<rank> neighbors;
  // Stores the "GhostExchangeDesc" elements of the process(partition)
  std::vector<GhostExchangeDesc> exchangeVector;

  globox::GlobalBox<int, int> gbox;

  // Stores the local cells and the ghost cells of this process(partition)
  std::vector<int> cells;
  // Stores all local cells which ar at the border of the local area of
  // this process (Stores the local index)
  std::vector<int> borderCells;
  // Stores for each cell in "borderCells" the ranks of their neighbourhood
  // Key is the local cell ID and value a set of ranks
  std::map<lidx, std::vector<rank>> borderCellsNeighbors;
  // Stores the received partition array
  // Can be dircetly used in Metis algorithm
  std::vector<rank> partition;
  // Stores the mapping of a global linearized index to a
  // local linearized index or an ghost linearized index.
  // Global index is key and local/ghost index is value.
  std::unordered_map<int, int> global_to_local;

  // After partition rebuild the local part of the structure
  void reinit(bool init = false);

  // Clears obsolete entries from "partition"
  void clear_unknown_cell_ownership();

  // Computes vector of vectors of cells which has to be send to neighbours
  std::vector<std::vector<int>>
  compute_send_list(std::vector<double>&& sendLoads,
                    const repart::Metric &m);

  MPI_Comm neighcomm;

  // Send message with neighbourhood of received cells in "sendCells"
  std::vector<std::vector<NeighSend>>
  sendNeighbourhood(const std::vector<std::vector<int>> &toSend);
  // Update partition array
  void updateReceivedNeighbourhood(
      const std::vector<std::vector<NeighSend>> &neighbourhood);
};
} // namespace grids
} // namespace generic_dd

#endif
