
#ifndef GENERIC_DD_GRIDS_GLOBOX_HPP_INCLUDED
#define GENERIC_DD_GRIDS_GLOBOX_HPP_INCLUDED

#include <array>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/iterator_range.hpp>
#include <cmath>

#include "cells.hpp" // max_range
#include "grid.hpp"  // node_grid, node_pos, box_l

namespace generic_dd {
namespace grids {
namespace globox {

template <typename GloBox>
struct NeighborIterator
    : public boost::iterator_facade<
          NeighborIterator<GloBox>,
          typename GloBox::index_type_1d,
          boost::random_access_traversal_tag,
          /* don't use reference type for return */
          typename GloBox::index_type_1d> {
private:
  using base_type = boost::iterator_facade<NeighborIterator<GloBox>,
                                           int,
                                           boost::random_access_traversal_tag,
                                           int>;
public:
  using value_type = typename base_type::value_type;
  using difference_type = typename base_type::difference_type;

  NeighborIterator() : cell(-1), idx(27), g(nullptr) {}
  NeighborIterator(const GloBox *g, value_type cell, int start)
      : cell(cell), idx(start), g(g) {}

private:
  friend class boost::iterator_core_access;
  //using difference_type = typename base_type::difference_type;

  value_type dereference() const { return g->neighbor(cell, idx); }

  bool equal(NeighborIterator const &other) const { return idx == other.idx; }

  void increment() { advance(1); }

  void decrement() { advance(-1); }

  void advance(difference_type n) { idx += n; }

  difference_type distance_to(NeighborIterator const &other) const {
    return other.idx - idx;
  }

  value_type cell;
  int idx;
  const GloBox *g;
};

inline const std::array<int, 3> &neigh_offset_3d(int i) {
  /*
   * Offsets of cell and process neighbors -- respects cell ordering
   * requirements of pargrid.hpp: 0 cell itself 1-13 hs neigh 14-26 fs neigh
   *
   * For iterating over all neighbors without 0, iterate over neigh_offset[i+1]
   * for i = 0..26
   */
  static constexpr std::array<std::array<int, 3>, 27> _neigh_offset_3d = {
      {{{0, 0, 0}},
       {{1, 0, 0}},
       {{-1, 1, 0}},
       {{0, 1, 0}},
       {{1, 1, 0}},
       {{-1, -1, 1}},
       {{0, -1, 1}},
       {{1, -1, 1}},
       {{-1, 0, 1}},
       {{0, 0, 1}},
       {{1, 0, 1}},
       {{-1, 1, 1}},
       {{0, 1, 1}},
       {{1, 1, 1}},
       // Full shell begin
       {{-1, -1, -1}},
       {{0, -1, -1}},
       {{1, -1, -1}},
       {{-1, 0, -1}},
       {{0, 0, -1}},
       {{1, 0, -1}},
       {{-1, 1, -1}},
       {{0, 1, -1}},
       {{1, 1, -1}},
       {{-1, -1, 0}},
       {{0, -1, 0}},
       {{1, -1, 0}},
       {{-1, 0, 0}}}};
  return _neigh_offset_3d[i];
}

template <typename T1, typename T2>
std::array<
  typename std::common_type<typename T1::value_type,
                            typename T2::value_type>::type,
  3> div_3(const T1& a, const T2& b) {
  return {{a[0] / b[0], a[1] / b[1], a[2] / b[2]}};
}

template <typename index1d, typename index3d = index1d> struct GlobalBox {
  typedef index1d index_type_1d;
  typedef index3d index_type_3d;
  typedef std::array<index_type_3d, 3> cell_index_type;
  typedef std::array<double, 3> position_type;

  cell_index_type m_cell_grid;
  cell_index_type m_cell_grid_corr;
  position_type m_cell_size;

  std::array<index_type_1d, 27> m_neigh_offset_1d;

  // Initialize with Espresso internals
  GlobalBox()
      : m_cell_grid({{static_cast<index_type_3d>(box_l[0] / max_range),
                      static_cast<index_type_3d>(box_l[1] / max_range),
                      static_cast<index_type_3d>(box_l[2] / max_range)}}),
        m_cell_size(div_3(box_l, m_cell_grid)) {
    m_cell_grid_corr[0] = linearize({m_cell_grid[0], 0, 0});
    m_cell_grid_corr[1] = linearize({0, m_cell_grid[1], 0});
    m_cell_grid_corr[2] = linearize({0, 0, m_cell_grid[2]});

    for (int i = 0; i < 27; ++i)
      m_neigh_offset_1d[i] = __linearize(neigh_offset_3d(i));
  }

  inline void apply_pbc(index_type_3d cell[3]) const noexcept {
    for (int d = 0; d < 3; d++) {
      cell[d] -= std::floor(cell[d] / static_cast<double>(m_cell_grid[d])) *
                 m_cell_grid[d];
    }
  }

  inline index_type_1d linearize(const index_type_3d *cell) const noexcept {
    return __linearize(cell);
  }

  inline index_type_1d linearize(const cell_index_type &cell) const noexcept {
    return __linearize(cell);
  }

  inline cell_index_type unlinearize(index_type_1d pos) const {
    cell_index_type idx;
    idx[2] = static_cast<index_type_3d>(pos % m_cell_grid[2]);
    pos /= m_cell_grid[2];
    idx[1] = static_cast<index_type_3d>(pos % m_cell_grid[1]);
    idx[0] = static_cast<index_type_3d>(pos / m_cell_grid[1]);
    return idx;
  }

  inline index_type_1d neighbor(index_type_1d index, int neigh) const {
    if (neigh < 0 || neigh >= 27)
      throw std::out_of_range("Neighbor index out of range.");
    auto idx = unlinearize(index);
    const auto &no = neigh_offset_3d(neigh);
    auto ni = index + m_neigh_offset_1d[neigh];

    for (int d = 0; d < 3; ++d) {
      // Can be out of bounds by at most 1, i.e. subtracting m_cell_grid once
      // suffices.
      if (idx[d] == 0 && no[d] < 0)
        ni += m_cell_grid_corr[d];
      else if (idx[d] == m_cell_grid[d] - 1 && no[d] > 0)
        ni -= m_cell_grid_corr[d];
    }
    return ni;
  }

  using NeighIt = NeighborIterator<GlobalBox>;
  boost::iterator_range<NeighIt> full_shell_neigh(index_type_1d index) {
    return {NeighIt(this, index, 0), NeighIt()};
  }

  boost::iterator_range<NeighIt> full_shell_neigh_without_center(index_type_1d index) {
    return {NeighIt(this, index, 1), NeighIt()};
  }

  inline index_type_1d cell_at_pos(const double pos[3]) const noexcept {
    cell_index_type cell;
    for (int d = 0; d < 3; d++) {
      cell[d] = static_cast<index_type_3d>(pos[d] / m_cell_size[d]);
    }

    return linearize(cell.data());
  }

  inline index_type_1d ncells() const noexcept {
    return static_cast<index_type_1d>(m_cell_grid[0]) * m_cell_grid[1] *
           m_cell_grid[2];
  }

  inline position_type cell_size() const noexcept {
    return m_cell_size;
  }

  inline cell_index_type grid_size() const noexcept { return m_cell_grid; }

  inline position_type midpoint(index_type_1d index) const {
    auto ii = unlinearize(index);
    position_type midpoint;

    for (int d = 0; d < 3; ++d)
      midpoint[d] = m_cell_size[d] * (ii[d] + .5);;
    return midpoint;
  }

private:
  // This is a template because we need it for index_type_3d and int.
  template <typename T>
  inline index_type_1d __linearize(const T *cell) const noexcept {
    return (static_cast<index_type_1d>(cell[0]) * m_cell_grid[1] + cell[1]) *
               m_cell_grid[2] +
           cell[2];
  }

  template <typename T>
  inline index_type_1d __linearize(const std::array<T, 3> &cell) const
      noexcept {
    return __linearize(cell.data());
  }
};

} // namespace globox
} // namespace grids
} // namespace generic_dd
#endif