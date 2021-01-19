/**
 * Copyright 2017-2020 Steffen Hirschmann
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

#include <unordered_map>
#include <vector>

#include "pargrid.hpp"

#include "ioptional.hpp"
#include "range.hpp"

namespace repa {
namespace util {

/** Class for storing local and ghost cell indices and the
 * mappings of global indices to local and ghost indices.
 *
 * Store is filled with cells vial push_back_local() and push_back_ghost().
 */
struct global_index_storage {

    /** Removes all elements
     */
    void clear()
    {
        _local_cells.clear();
        _ghost_cells.clear();
        _inverse_map.clear();
    }

    /** Registers a new local cell index
     *
     * @param g global index of the new local cell
     */
    void push_back_local(global_cell_index_type g)
    {
        assert(!holds_global_index(g));
        _inverse_map.emplace(
            g, local_cell_index_type{
                   _local_cells.size()}); // Enumerate the cells from 0 to n
        _local_cells.push_back(g);
    }

    /** Alternative to "push_back_local" if all local cells are given
     * in an array.
     */
    void bulk_register_local_cells(std::vector<global_cell_index_type> &&cs)
    {
        assert(_local_cells.empty());
        _local_cells = std::move(cs);
        size_t i = 0;
        for (const auto &g : _local_cells) {
            _inverse_map.emplace(g, local_cell_index_type{i++});
        }
    }

    /** Registers a new ghost cell index
     *
     * @param g global index of the new ghost cell
     */
    void push_back_ghost(global_cell_index_type g)
    {
        assert(!holds_global_index(g));
        _inverse_map.emplace(
            g, ghost_cell_index_type{
                   _ghost_cells.size()}); // Enumerate the cells from 0 to m
        _ghost_cells.push_back(g);
    }

    /** Returns the associated global index of a local cell index
     */
    global_cell_index_type as_global_index(local_cell_index_type index) const
    {
        assert(index >= 0 && static_cast<size_t>(index) < _local_cells.size());
        return _local_cells[index];
    }

    /** Returns the associated global index of a ghost cell index.
     */
    global_cell_index_type as_global_index(ghost_cell_index_type index) const
    {
        assert(index >= 0 && static_cast<size_t>(index) < _ghost_cells.size());
        return _ghost_cells[index];
    }

    /** Returns the associated global index of a local or a ghost cell index.
     */
    global_cell_index_type
    as_global_index(local_or_ghost_cell_index_type index) const
    {
        auto apply_operator_at
            = [this](const auto &v) { return this->as_global_index(v); };
        return index.fmap(apply_operator_at);
    }

    /** Returns the local or ghost index of a global cell index.
     */
    local_or_ghost_cell_index_type
    as_local_or_ghost_index(global_cell_index_type g) const
    {
        return _inverse_map.at(g);
    }

    /** Returns true if this store holds a global cell index as local or
     * ghost index.
     */
    bool holds_global_index(global_cell_index_type g) const
    {
        return _inverse_map.find(g) != std::end(_inverse_map);
    }

    /** Returns a local cell index of a global cell index. If the global cell
     * index is not associated with a local index, returns an empty value.
     */
    ioptional<local_cell_index_type>
    as_local_index(global_cell_index_type g) const
    {
        auto lgidx = as_local_or_ghost_index(g);
        if (lgidx.is<local_cell_index_type>())
            return lgidx.as<local_cell_index_type>();
        else
            return {};
    }

    /** Returns a ghost cell index of a global cell index. If the global cell
     * index is not associated with a ghost index, returns an empty value.
     */
    ioptional<ghost_cell_index_type>
    as_ghost_index(global_cell_index_type g) const
    {
        auto lgidx = as_local_or_ghost_index(g);
        if (lgidx.is<ghost_cell_index_type>())
            return lgidx.as<ghost_cell_index_type>();
        else
            return {};
    }

    /** Returns a range of all local cells
     */
    auto local_cells() const
    {
        return util::range(local_cell_index_type{_local_cells.size()});
    }

    /** Returns a range of all ghost cells
     */
    auto ghost_cells() const
    {
        return util::range(ghost_cell_index_type{_ghost_cells.size()});
    }

    /** Returns a lambda that resolves global to local positions
     * while ensuring that the indices are purely local ones.
     */
    auto global_to_local_transformer() const
    {
        return [this](const global_cell_index_type gloidx) {
            return as_local_index(gloidx).ensure_value();
        };
    }

    /** Returns a lambda that resolves global to ghost positions
     * while ensuring that the indices are purely ghost ones.
     */
    auto global_to_ghost_transformer() const
    {
        return [this](const global_cell_index_type gloidx) {
            return as_ghost_index(gloidx).ensure_value();
        };
    }

private:
    std::vector<global_cell_index_type> _local_cells;
    std::vector<global_cell_index_type> _ghost_cells;

    std::unordered_map<global_cell_index_type, local_or_ghost_cell_index_type>
        _inverse_map;
};

} // namespace util
} // namespace repa