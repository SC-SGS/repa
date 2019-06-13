
#pragma once

//#ifdef HAVE_P4EST

#include "util/p4est_deleter.hpp"

#include <array>
#include <vector>
#include <memory>
#include "../pargrid.hpp"



namespace repa {
namespace grids {

namespace impl {

enum class CellType { inner = 0, boundary = 1, ghost = 2 };
struct LocalShell {
  int idx; // a unique index within all cells (as used by p8est for locals)
  int rank; // the rank of this cell (equals this_node for locals)
  CellType shell; // shell information (0: inner local cell, 1: boundary local cell, 2: ghost cell)
  int boundary; // Bit mask storing boundary info. MSB ... z_r,z_l,y_r,y_l,x_r,x_l LSB
                // Cells with shell-type 0 or those located within the domain are always 0
                // Cells with shell-type 1 store information about which face is a boundary
                // Cells with shell-type 2 are set if the are in the periodic halo
  std::array<int, 26> neighbor; // unique index of the fullshell neighborhood cells (as in p8est); only 26 because cell itself is not included.
  Vec3i coord; // cartesian coordinates of the cell

  LocalShell(int idx, int rank, CellType shell, int boundary, int x, int y,
             int z)
      : idx(idx), rank(rank), shell(shell), boundary(boundary), coord({{x, y, z}}) {
    neighbor.fill(-1);
  }
};

struct RepartState {
  const boost::mpi::communicator& comm_cart;
  bool after_repart;
  std::vector<p4est_locidx_t> nquads_per_proc;
  std::function<void()> exchange_start_callback;

  RepartState(const boost::mpi::communicator& comm_cart);

  inline void reset();
  inline void inc_nquads(rank proc);
  inline void allreduce();
};

}

struct P4estGrid : public ParallelLCGrid {
  P4estGrid(const boost::mpi::communicator& comm, Vec3d box_size, double min_cell_size);
  lidx n_local_cells() override;
  gidx n_ghost_cells() override;
  nidx n_neighbors() override;
  rank neighbor_rank(nidx i) override;
  std::array<double, 3> cell_size() override;
  Vec3i grid_size() override;
  lgidx cell_neighbor_index(lidx cellidx, int neigh) override;
  std::vector<GhostExchangeDesc> get_boundary_info() override;
  lidx position_to_cell_index(double pos[3]) override;
  rank position_to_rank(double pos[3]) override;
  nidx position_to_neighidx(double pos[3]) override;
  bool repartition(const repart::Metric& m, std::function<void()> exchange_start_callback) override;

private:
  int m_grid_level;

  // Number of grid cells in total and per tree.
  Vec3i m_grid_size, m_brick_size;
  // Cell size (box_l / m_grid_size)
  std::array<double, 3> m_cell_size, m_inv_cell_size;

  // p4est data structures
  std::unique_ptr<p8est_connectivity_t> m_p8est_conn;
  std::unique_ptr<p8est_t> m_p8est;
  int m_num_local_cells, m_num_ghost_cells;

  // helper data structures
  std::vector<int> m_node_first_cell_idx;
  std::vector<impl::LocalShell> m_p8est_shell;

  // comm data structures
  std::vector<GhostExchangeDesc> m_exdescs;
  std::vector<int> m_neighranks;

  void set_optimal_cellsize();
  void create_grid();
  void prepare_communication();
  // Reinitialized the grid (instantiation or after repartitioning)
  void reinitialize();

  impl::RepartState m_repartstate;
};
}
}

//#endif // HAVE_P4EST