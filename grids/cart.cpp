
#include "cart.hpp"
#include "communication.hpp" // comm_cart
#include "grid.hpp" // node_grid, node_pos, box_l
#include "domain_decomposition.hpp"

#include <algorithm>

#define UNUSED(x) ((void) (x))

namespace generic_dd {
namespace grids {

namespace impl {

template <typename T>
void push_back_unique(std::vector<T>& v, T val)
{
  if (std::find(std::begin(v), std::end(v), val) == std::end(v))
    v.push_back(val);
}

// Offsets of cell and process neighbors -- respects cell ordering requirements
// of pargrid.hpp:
// 0 cell itself
// 1-13 hs neigh
// 14-26 fs neigh
static std::vector<std::array<int, 3>> cart_neigh_offset = {
  { 0, 0, 0},
  { 1, 0, 0},
  {-1, 1, 0},
  { 0, 1, 0},
  { 1, 1, 0},
  {-1,-1, 1},
  { 0,-1, 1},
  { 1,-1, 1},
  {-1, 0, 1},
  { 0, 0, 1},
  { 1, 0, 1},
  {-1, 1, 1},
  { 0, 1, 1},
  { 1, 1, 1},
  // Full shell begin
  {-1,-1,-1},
  { 0,-1,-1},
  { 1,-1,-1},
  {-1, 0,-1},
  { 0, 0,-1},
  { 1, 0,-1},
  {-1, 1,-1},
  { 0, 1,-1},
  { 1, 1,-1},
  {-1,-1, 0},
  { 0,-1, 0},
  { 1,-1, 0},
  {-1, 0, 0}
};

static int linearize(const std::array<int, 3>& c, const std::array<int, 3>& grid)
{
  return (c[0] * grid[1] + c[1]) * grid[2] + c[2];
}

static std::array<int, 3> unlinearize(int cidx, const std::array<int, 3>& grid)
{
  return {{ (cidx / grid[2]) / grid[1],
            (cidx / grid[2]) % grid[1],
            cidx % grid[2] }};
}

static std::array<int, 3> pos_add_folded(const std::array<int, 3>& pos, const std::array<int, 3>& offset, const std::array<int, 3>& grid)
{
  std::array<int, 3> res;

  for (int i = 0; i < 3; ++i) {
    res[i] = pos[i] + offset[i];
    if (res[i] < 0 || res[i] >= grid[i])
      res[i] -= (res[i] / grid[i]) * grid[i];
  }
  return res;
}


static std::pair<std::array<int, 3>, std::array<int, 3>>
determine_send_receive_bounds(const std::array<int, 3>& offset, int receive, const std::array<int, 3>& grid)
{
  std::array<int, 3> lc, hc;

  for (int i = 0; i < 3; ++i) {
    lc[i] = offset[i] <= 0? 1: grid[i];
    hc[i] = offset[i] < 0? 1: grid[i];

    // The receive area is actually in the ghost layer
    // so shift the corresponding indices.
    if (receive) {
      if (offset[i] > 0)
        lc[i] = hc[i] = grid[i] + 1;
      else if (offset[i] < 0)
        lc[i] = hc[i] = 0;
    }
  }
  return std::make_pair(lc, hc);
}

}

bool CartGrid::is_ghost_cell(const std::array<int, 3>& c)
{
  return c[0] == 0 || c[0] == m_ghost_grid_size[0] - 1 ||
         c[1] == 0 || c[1] == m_ghost_grid_size[1] - 1 ||
         c[2] == 0 || c[2] == m_ghost_grid_size[2] - 1;
}

rank CartGrid::proc_offset_to_rank(const std::array<int, 3> &offset)
{
  auto neighpos = impl::pos_add_folded(m_procgrid_pos, offset, m_procgrid);
  int rank;
  MPI_Cart_rank(comm_cart, neighpos.data(), &rank);
  return rank;
}

void CartGrid::fill_neighranks()
{
  m_neighranks.clear();

  for (const auto& offset: impl::cart_neigh_offset) {
    // Push back unique neighbor ranks into m_neighbors
    impl::push_back_unique(m_neighranks, proc_offset_to_rank(offset));
  }
}

void CartGrid::create_index_permutations()
{
  int ncells = n_local_cells() + n_ghost_cells();
  m_to_pargrid_order.resize(ncells);
  m_from_pargrid_order.resize(ncells);

  int localidx = 0, ghostidx = n_local_cells();
  for (int i = 0; i < ncells; ++i) {
    auto c = impl::unlinearize(i, m_ghost_grid_size);
    if (is_ghost_cell(c)) {
      m_from_pargrid_order[ghostidx] = i;
      m_to_pargrid_order[i] = ghostidx++;
    } else {
      m_from_pargrid_order[localidx] = i;
      m_to_pargrid_order[i] = localidx++;
    }
  }
}

void CartGrid::create_grid()
{
  for (int i = 0; i < 3; ++i) {
    // Copy infos from grid.hpp (dependant on comm_cart)
    m_procgrid[i] = node_grid[i];
    m_procgrid_pos[i] = node_pos[i];

    // Local box info
    m_localbox[i] = box_l[i] / m_procgrid[i];
    m_lowerleft[i] = m_localbox[i] * m_procgrid_pos[i];

    // Grid and cell size
    if (max_range > ROUND_ERROR_PREC * box_l[0])
      m_grid_size[i] = static_cast<int>(m_localbox[i] / max_range);
    else
      m_grid_size[i] = 1;
    m_ghost_grid_size[i] = m_grid_size[i] + 2;

    m_cell_size[i] = m_localbox[i] / m_grid_size[i];
    m_inv_cell_size[i] = 1.0 / m_cell_size[i];
  }
}


void CartGrid::fill_comm_cell_lists(std::vector<int>& v, const std::array<int, 3>& lc, const std::array<int, 3>& hc)
{
  std::array<int, 3> c;
  for (c[0] = lc[0]; c[0] <= hc[0]; c[0]++) {
    for (c[1] = lc[1]; c[1] <= hc[1]; c[1]++) {
      for (c[2] = lc[2]; c[2] <= hc[2]; c[2]++) {
        v.push_back(linearize(c));
      }
    }
  }
}

void CartGrid::prepare_communication()
{
  m_exdescs.clear();
  m_exdescs.resize(n_neighbors());

  // The loop below is not guaranteed to get to the node itself (if it is
  // neighbored in every direction by other processes). Therefore, set this
  // communication destination fix beforehand.
  m_exdescs[neighbor_idx(this_node)].dest = this_node;

  for (auto o = 1; o < impl::cart_neigh_offset.size(); ++o) {
    std::array<int, 3> lc, hc;
    const auto& offset = impl::cart_neigh_offset[o];
    std::array<int, 3> opposite = {{ -offset[0], -offset[1], -offset[2] }};

    // Send
    {
      auto neigh = proc_offset_to_rank(offset);
      auto i = neighbor_idx(neigh);
      m_exdescs[i].dest = neigh;
      std::tie(lc, hc) = impl::determine_send_receive_bounds(offset, 0, m_grid_size);
      fill_comm_cell_lists(m_exdescs[i].send, lc, hc);
    }

    // Receive in opposite direction. Otherwise send and receive order on
    // the processes won't match.
    {
      auto neigh = proc_offset_to_rank(opposite);
      auto i = neighbor_idx(neigh);
      m_exdescs[i].dest = neigh;
      std::tie(lc, hc) = impl::determine_send_receive_bounds(opposite, 1, m_grid_size);
      fill_comm_cell_lists(m_exdescs[i].recv, lc, hc);
    }
  }
}

CartGrid::CartGrid()
{
  create_grid();
  create_index_permutations();
  fill_neighranks();
  prepare_communication();
}

lidx CartGrid::n_local_cells()
{
  return m_grid_size[0] * m_grid_size[1] * m_grid_size[2];
}

gidx CartGrid::n_ghost_cells()
{
  int ggs = m_ghost_grid_size[0] * m_ghost_grid_size[1] * m_ghost_grid_size[2];
  return ggs - n_local_cells();
}

nidx CartGrid::n_neighbors()
{
  return m_neighranks.size();
}

rank CartGrid::neighbor_rank(nidx i)
{
  return m_neighranks[i];
}

lgidx CartGrid::cell_neighbor_index(lidx cellidx, int neigh)
{
  auto c = unlinearize(cellidx);
  auto nc = impl::pos_add_folded(c, impl::cart_neigh_offset[neigh],
                                 m_ghost_grid_size);
  return linearize(nc);
}

lgidx CartGrid::linearize(std::array<int, 3> c)
{
  auto idx = impl::linearize(c, m_ghost_grid_size);
  return m_to_pargrid_order[idx];
}

std::array<int, 3> CartGrid::unlinearize(lgidx cidx)
{
  auto idx = m_from_pargrid_order[cidx];
  return impl::unlinearize(idx, m_ghost_grid_size);
}

std::vector<GhostExchangeDesc> CartGrid::get_boundary_info()
{
  return m_exdescs;
}

lidx CartGrid::position_to_cell_index(double pos[3])
{
  std::array<int, 3> c;

  for (int i = 0; i < 3; ++i) {
    // Transform to process local coordinates
    double tpos = pos[i] - m_lowerleft[i];
    if (tpos < 0.0 || tpos >= m_localbox[i])
      throw std::domain_error("Particle not in local box");

    c[i] = tpos * m_inv_cell_size[i] + 1; // +1 to skip the ghost cells
  }
  return linearize(c);
}

rank CartGrid::position_to_rank(double pos[3])
{
  std::array<int, 3> proc;
  for (int i = 0; i < 3; ++i)
    proc[i] = static_cast<int>(pos[i] * m_inv_cell_size[i]) / m_grid_size[i];

  int rank;
  MPI_Cart_rank(comm_cart, proc.data(), &rank);
  return rank;
}

nidx CartGrid::neighbor_idx(rank r)
{
  // Search this rank in the local neighbor list and return its index
  // Use std::find here as 1) m_neighranks might not be sorted and 2) it has
  // at most 26 entries, so sequential search might not hurt that much.
  auto it = std::find(std::begin(m_neighranks), std::end(m_neighranks), r);
  if (*it != r)
    throw std::runtime_error("Rank not a neighbor.");

  return std::distance(std::begin(m_neighranks), it);
}

nidx CartGrid::position_to_neighidx(double pos[3])
{
  // Determine the neighbor rank for locally known cell
  // Using position_to_rank here as it is the simpler code. Could also
  // search the neighboring cells of the cell where pos lies in.
  auto rank = position_to_rank(pos);
  return neighbor_idx(rank);
}

std::array<double, 3> CartGrid::cell_size()
{
  return m_cell_size;
}

std::array<int, 3> CartGrid::grid_size()
{
  return m_grid_size;
}

bool CartGrid::repartition(const repart::Metric& m,
                           std::function<void()> exchange_start_callback)
{
  UNUSED(m);
  UNUSED(exchange_start_callback);
  return false;
}

}
}
