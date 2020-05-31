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

#pragma once

#include <array>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/iterator_range.hpp>
#include <cmath>

#include "_compat.hpp"
#include "util/linearize.hpp"
#include "util/neighbor_offsets.hpp"
#include "util/range.hpp"

namespace repa {
namespace grids {
namespace globox {

/** Random access iterator class for iterating over neighborhoods of cells.
 */
template <typename GloBox>
struct NeighborIterator
    : public boost::iterator_facade<NeighborIterator<GloBox>,
                                    typename GloBox::index_type_1d,
                                    boost::random_access_traversal_tag,
                                    /* don't use reference type for return */
                                    typename GloBox::index_type_1d> {
private:
    using base_type = boost::iterator_facade<NeighborIterator<GloBox>,
                                             typename GloBox::index_type_1d,
                                             boost::random_access_traversal_tag,
                                             typename GloBox::index_type_1d>;

public:
    using value_type = typename base_type::value_type;
    using difference_type = typename base_type::difference_type;

    NeighborIterator() = delete;

    NeighborIterator(const GloBox *g, value_type cell, int start)
        : cell(cell), idx(start), g(g)
    {
    }

private:
    friend class boost::iterator_core_access;
    // using difference_type = typename base_type::difference_type;

    value_type dereference() const
    {
        return g->neighbor(cell, idx);
    }

    bool equal(NeighborIterator const &other) const
    {
        return idx == other.idx;
    }

    void increment()
    {
        advance(1);
    }

    void decrement()
    {
        advance(-1);
    }

    void advance(difference_type n)
    {
        idx += n;
    }

    difference_type distance_to(NeighborIterator const &other) const
    {
        return other.idx - idx;
    }

    value_type cell;
    difference_type idx;
    const GloBox *g;
};

/** Global cell ordering.
 *
 */
template <typename index1d, typename index3d = index1d>
struct GlobalBox {
    typedef index1d index_type_1d;
    typedef index3d index_type_3d;
    typedef Vec3<index_type_3d> cell_index_type;
    typedef Vec3d position_type;

    cell_index_type m_cell_grid;
    cell_index_type m_cell_grid_corr;
    position_type m_cell_size;

    std::array<index_type_1d, 27> m_neigh_offset_1d;

    /** Initialization with a regular grid
     * @param box_size box size
     * @param min_cell_size minimum cell size
     */
    GlobalBox(Vec3d box_size, double min_cell_size)
        : m_cell_grid(static_cast<index_type_3d>(box_size[0] / min_cell_size),
                      static_cast<index_type_3d>(box_size[1] / min_cell_size),
                      static_cast<index_type_3d>(box_size[2] / min_cell_size)),
          m_cell_size(box_size[0] / m_cell_grid[0],
                      box_size[1] / m_cell_grid[1],
                      box_size[2] / m_cell_grid[2])
    {
        constexpr auto zero = index_type_3d{0};
        m_cell_grid_corr[0]
            = linearize(cell_index_type{m_cell_grid[0], zero, zero});
        m_cell_grid_corr[1]
            = linearize(cell_index_type{zero, m_cell_grid[1], zero});
        m_cell_grid_corr[2]
            = linearize(cell_index_type{zero, zero, m_cell_grid[2]});

        std::transform(std::begin(util::NeighborOffsets3D::raw),
                       std::end(util::NeighborOffsets3D::raw),
                       std::begin(m_neigh_offset_1d), [this](const auto &cell) {
                           return this->linearize(cell);
                           //     ^^^^^^ Note: gcc-5.x compatibility
                       });
    }

    /** Returns the index of a neighboring cell
     * @param index index of the center cell
     * @param neigh index [0, 27) of the neighbor
     */
    inline index_type_1d neighbor(index_type_1d index, fs_neighidx neigh) const
    {
        auto idx = unlinearize(index);
        const auto &no = util::NeighborOffsets3D::raw[neigh];
        auto ni = index + m_neigh_offset_1d[neigh];

        for (int d = 0; d < 3; ++d) {
            // Can be out of bounds by at most 1, i.e. subtracting m_cell_grid
            // once suffices.
            if (idx[d] == 0 && no[d] < 0)
                ni += index_type_1d{m_cell_grid_corr[d]};
            else if (idx[d] == m_cell_grid[d] - 1 && no[d] > 0)
                ni -= index_type_1d{m_cell_grid_corr[d]};
        }
        return ni;
    }

    /** Returns a range object that allows iterating over the full
     * shell neighborhood (including the center cell "index" itself)
     */
    using NeighIt = NeighborIterator<GlobalBox>;
    boost::iterator_range<NeighIt> full_shell_neigh(index_type_1d index) const
    {
        return {NeighIt(this, index, 0), NeighIt(this, index, 27)};
    }

    /** Returns a range object that allows iterating over the full
     * shell neighborhood without the center cell ("index").
     */
    boost::iterator_range<NeighIt>
    full_shell_neigh_without_center(index_type_1d index) const
    {
        return {NeighIt(this, index, 1), NeighIt(this, index, 27)};
    }

    /** Returns the index of the cell at position "pos".
     * @param pos coordinates of the position in [0.0, box_size[i])
     */
    inline index_type_1d cell_at_pos(Vec3d pos) const noexcept
    {
        cell_index_type cell;
        for (int d = 0; d < 3; d++) {
            cell[d] = static_cast<index_type_3d>(pos[d] / m_cell_size[d]);
        }

        return linearize(cell);
    }

    /** Returns an iterator range to all global cell indices
     */
    inline auto global_cells() const
    {
        return util::range(ncells());
    }

    inline bool is_valid_global_index(index_type_1d i) const
    {
        return i >= 0 && i < ncells();
    }

    /** Returns the resulting cell size, greater or equal to min_cell_size.
     */
    inline position_type cell_size() const noexcept
    {
        return m_cell_size;
    }

    /** Returns the number of cells in each dimension
     */
    inline cell_index_type grid_size() const noexcept
    {
        return m_cell_grid;
    }

    /** Returns the midpoint of a cell
     */
    inline position_type midpoint(index_type_1d index) const
    {
        auto ii = unlinearize(index);
        position_type midpoint;

        for (int d = 0; d < 3; ++d)
            midpoint[d] = m_cell_size[d] * (ii[d] + .5);
        return midpoint;
    }

private:
    // This is a template because we need it for index_type_3d and int.
    template <typename T>
    inline index_type_1d linearize(const Vec3<T> &cell) const noexcept
    {
        return index_type_1d{util::linearize(cell, m_cell_grid)};
    }

    inline cell_index_type unlinearize(index_type_1d pos) const
    {
        return util::unlinearize(pos, m_cell_grid);
    }

    /** Returns the number of cells
     */
    inline index_type_1d ncells() const noexcept
    {
        return index_type_1d{m_cell_grid[0] * m_cell_grid[1] * m_cell_grid[2]};
    }

}; // namespace globox

} // namespace globox
} // namespace grids
} // namespace repa