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

#include "grid_types.hpp"

namespace repa {

namespace {

static const std::unordered_map<std::string, GridType> grid_type_registry
    = {{"p4est", GridType::P4EST},
       {"cart", GridType::CART},
       {"graph", GridType::GRAPH},
       {"diff", GridType::DIFF},
       {"hybrid_gp_diff", GridType::HYB_GP_DIFF},
       {"kd_tree", GridType::KD_TREE},
       {"gridbased", GridType::GRIDBASED}};

#ifdef HAVE_KDPART
#define KDPART_AVAIL true
#else
#define KDPART_AVAIL false
#endif

#ifdef HAVE_P4EST
#define P4EST_AVAIL true
#else
#define P4EST_AVAIL false
#endif

#ifdef HAVE_METIS
#define METIS_AVAIL true
#else
#define METIS_AVAIL false
#endif

#ifdef HAVE_TETRA
#define TETRA_AVAIL true
#else
#define TETRA_AVAIL false
#endif

// Awaiting C++20 which will finally have designated initializers -.-
static const std::unordered_map<GridType, bool> grid_type_availability
    = {{GridType::P4EST, P4EST_AVAIL},
       {GridType::CART, true},
       {GridType::GRAPH, METIS_AVAIL},
       {GridType::DIFF, true},
       {GridType::HYB_GP_DIFF, METIS_AVAIL},
       {GridType::KD_TREE, KDPART_AVAIL},
       {GridType::GRIDBASED, true}};

} // namespace

inline GridType parse_grid_type(const std::string &desc)
{
    try {
        return grid_type_registry.at(desc);
    }
    catch (const std::out_of_range &) {
        throw UnknownGridTypeError(desc);
    }
}

std::string grid_type_to_string(GridType gt)
{
    for (const auto &p : grid_type_registry) {
        if (gt == p.second)
            return p.first;
    }

    throw UnknownGridTypeError();
}

bool has_grid_type(GridType gt)
{
    try {
        return grid_type_availability.at(gt);
    }
    catch (const std::out_of_range &) {
        throw UnknownGridTypeError();
    }
}

std::set<GridType> supported_grid_types()
{
    std::set<GridType> gts;
    for (const auto &p : grid_type_availability) {
        if (p.second)
            gts.insert(p.first);
    }
    return gts;
}

} // namespace repa