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

#include <set>
#include <string>

namespace repa {

/** Enum of supported grid types.
 * Caution: Might not all be compiled in.
 */
enum class GridType {
    NONE,
    P4EST,
    CART,
    GRAPH,
    DIFF,
    KD_TREE,
    HYB_GP_DIFF,
    GRIDBASED
};

struct UnknownGridTypeError {
    UnknownGridTypeError() : w(std::string("Unknown grid type."))
    {
    }
    UnknownGridTypeError(std::string s)
        : w(std::string("Unknown grid type: `") + s + std::string("'"))
    {
    }
    virtual const char *what() const
    {
        return w.c_str();
    }

private:
    std::string w;
};

/** Returns the GridType associated with a descriptive string for the grid
 * type.
 */
GridType parse_grid_type(const std::string &desc);

/** Inverse of 'parse_grid_type'
 */
std::string grid_type_to_string(GridType gt);

/** Returns true if support for a certain grid type is compiled in.
 */
bool has_grid_type(GridType gt);

/** Returns a set of all supported grid types.
 */
std::set<GridType> supported_grid_types();

} // namespace repa
