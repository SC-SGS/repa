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
#include "linearize.hpp"
#include "range.hpp"
#include "vec_arith.hpp"


namespace repa {
namespace util {

#define _REPA_NDIR(d) ((d + 1) % 3)
#define _REPA_NNDIR(d) ((d + 2) % 3)

struct box_global_index_storage {
    typedef int local_box_index_type;

    /** Removes all elements
     */
    void clear()
    {
        _ghost_cells.clear();
        _inverse_ghost_map.clear();
    }

    void init(Vec3i global_grid_size,
              Vec3i local_grid_size,
              Vec3i local_lower_left_idx3d)
    {
        using namespace util::vector_arithmetic;

        _global_grid_size = global_grid_size;
        _local_lower_left_idx3d = local_lower_left_idx3d;
        _local_grid_size = local_grid_size;

        // Reserve full halo
        _ghost_cells.reserve(product(_local_grid_size + 2)
                             - product(_local_grid_size));

        auto ghost_cell_no = ghost_cell_index_type{0};
        // Assemble ghost->global and global->ghost mappings -- iterate over
        // full halo and find canonical representation of cell
        // Careful: this gets ugly.
        // Halo are cells that have at least one
        // coordinate == -1 or == local_grid_size[d].
        for (int d = 0; d < 3; ++d) {
            Vec3i local_ghost_idx3d;
            // The outer most loop has only two iterations
            // left (-1) and right (_local_grid_size[d]) side of the volume
            // spanned by _local_grid_size
            for (local_ghost_idx3d[d] = -1;
                 local_ghost_idx3d[d] <= _local_grid_size[d];
                 local_ghost_idx3d[d] += _local_grid_size[d] + 1) {

                // Iterate over 2d plane in codimension(d) with coordinate
                // ghost_idx3d[d] fixed.
                for (local_ghost_idx3d[_REPA_NDIR(d)] = -1;
                     local_ghost_idx3d[_REPA_NDIR(d)]
                     <= _local_grid_size[_REPA_NDIR(d)];
                     local_ghost_idx3d[_REPA_NDIR(d)]++) {

                    for (local_ghost_idx3d[_REPA_NNDIR(d)] = -1;
                         local_ghost_idx3d[_REPA_NNDIR(d)]
                         <= _local_grid_size[_REPA_NNDIR(d)];
                         local_ghost_idx3d[_REPA_NNDIR(d)]++) {

                        //const Vec3i gg
                        //    = _local_lower_left_idx3d + local_ghost_idx3d;
                        // std::cout << "At cell {" << gg[0] << "," << gg[1] <<
                        // ","
                        //          << gg[2] << "} -> r =";

                        const auto canon_idx3d = get_canonical_vector(
                            _local_lower_left_idx3d + local_ghost_idx3d);

                        // Sort out local cells
                        if (all(canon_idx3d >= _local_lower_left_idx3d)
                            && all(canon_idx3d < (_local_lower_left_idx3d
                                                  + _local_grid_size)))
                            continue;

                        const auto canon_idx
                            = get_global_cell_index(canon_idx3d);

                        // std::cout << " " << canon_idx.value() << std::endl;

                        if (_inverse_ghost_map.find(canon_idx)
                            != std::end(_inverse_ghost_map))
                            continue;

                        // std::cout << " > Inserting: " << canon_idx.value() <<
                        // ""
                        //          << std::endl;

                        _inverse_ghost_map.emplace(canon_idx, ghost_cell_no);
                        _ghost_cells.push_back(canon_idx);
                        ghost_cell_no++;
                    }
                }
            }
        }

        // Free difference between full halo and actual number of ghost cells
        _ghost_cells.shrink_to_fit();
    }

    global_cell_index_type as_global_index(local_cell_index_type index) const
    {
        using namespace util::vector_arithmetic;
        const Vec3i local_idx3d = util::unlinearize(index, _local_grid_size);
        return get_global_cell_index(local_idx3d + _local_lower_left_idx3d);
    }

    global_cell_index_type as_global_index(ghost_cell_index_type index) const
    {
        assert(index >= 0 && static_cast<size_t>(index) < _ghost_cells.size());
        return _ghost_cells[index];
    }

    global_cell_index_type
    as_global_index(local_or_ghost_cell_index_type index) const
    {
        auto apply_operator_at
            = [this](const auto &v) { return this->as_global_index(v); };
        return index.fmap(apply_operator_at);
    }

    local_or_ghost_cell_index_type
    as_local_or_ghost_index(global_cell_index_type g) const
    {
        using namespace util::vector_arithmetic;
        const Vec3i global_idx3d = util::unlinearize(g, _global_grid_size);

        if (all(global_idx3d >= _local_lower_left_idx3d)
            && all(global_idx3d
                   < (_local_lower_left_idx3d + _local_grid_size))) {
            // Local cell
            const Vec3i local_idx3d = global_idx3d - _local_lower_left_idx3d;
            return local_cell_index_type{
                util::linearize(local_idx3d, _local_grid_size)};
        }
        else {
            // Ghost cell
            assert(_inverse_ghost_map.find(g) != std::end(_inverse_ghost_map));
            return _inverse_ghost_map.at(g);
        }
    }

    // bool holds_global_index(global_cell_index_type g) const
    // {
    //     using namespace util::vector_arithmetic;
    //     const Vec3i global_idx3d = util::unlinearize(g, _global_grid_size);

    //     return (all(global_idx3d >= _local_lower_left_idx3d)
    //         && all(global_idx3d
    //                < (_local_lower_left_idx3d + _local_grid_size))) ||
    //         return _inverse_ghost_map.find(g) !=
    //         std::end(_inverse_ghost_map);
    // }

    ioptional<local_cell_index_type>
    as_local_index(global_cell_index_type g) const
    {
        auto lgidx = as_local_or_ghost_index(g);
        if (lgidx.is<local_cell_index_type>())
            return lgidx.as<local_cell_index_type>();
        else
            return {};
    }

    ioptional<ghost_cell_index_type>
    as_ghost_index(global_cell_index_type g) const
    {
        auto lgidx = as_local_or_ghost_index(g);
        if (lgidx.is<ghost_cell_index_type>())
            return lgidx.as<ghost_cell_index_type>();
        else
            return {};
    }

    ioptional<local_cell_index_type>
    as_local_index(const Vec3i &global_idx3d) const
    {
        return as_local_index(get_canonical_representant(global_idx3d));
    }

    ioptional<ghost_cell_index_type>
    as_ghost_index(const Vec3i &global_idx3d) const
    {
        return as_ghost_index(get_canonical_representant(global_idx3d));
    }

    auto local_cells() const
    {
        using namespace util::vector_arithmetic;
        return util::range(local_cell_index_type{product(_local_grid_size)});
    }

    auto ghost_cells() const
    {
        return util::range(ghost_cell_index_type{_ghost_cells.size()});
    }

    // /** Returns a lambda that resolves global to local positions
    //  * while ensuring that the indices are purely local ones.
    //  */
    // auto global_to_local_transformer() const
    // {
    //     return [this](const global_cell_index_type gloidx) {
    //         return as_local_index(gloidx).ensure_value();
    //     };
    // }

    // /** Returns a lambda that resolves global to ghost positions
    //  * while ensuring that the indices are purely ghost ones.
    //  */
    // auto global_to_ghost_transformer() const
    // {
    //     return [this](const global_cell_index_type gloidx) {
    //         return as_ghost_index(gloidx).ensure_value();
    //     };
    // }

    global_cell_index_type
    get_canonical_representant(const Vec3i &global_idx3d) const
    {
        return get_global_cell_index(get_canonical_vector(global_idx3d));
    }

private:
    Vec3i _global_grid_size;
    Vec3i _local_grid_size;
    Vec3i _local_lower_left_idx3d;

    std::vector<global_cell_index_type> _ghost_cells;

    std::unordered_map<global_cell_index_type, ghost_cell_index_type>
        _inverse_ghost_map;

    Vec3i get_canonical_vector(const Vec3i &global_idx3d) const
    {
        using namespace util::vector_arithmetic;
        return global_idx3d % _global_grid_size;
    }

    global_cell_index_type
    get_global_cell_index(const Vec3i &global_idx3d) const
    {
        return global_cell_index_type{
            util::linearize(global_idx3d, _global_grid_size)};
    }
};

} // namespace util
} // namespace repa