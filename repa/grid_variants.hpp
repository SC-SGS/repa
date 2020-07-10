/**
 * Copyright 2017-2020 The repa authors
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
#include <stdexcept>
#include <string>

#include "pargrid.hpp"

namespace repa {

/** Interface to query and set variants.
 */
struct VariantSetter {
    virtual std::set<std::string> get_supported_variants() const = 0;
    virtual void set_variant(const std::string &) = 0;
};

struct UnknownVariantError : public std::runtime_error {
    UnknownVariantError() : std::runtime_error("Unknown grid variant.")
    {
    }
};

/** Returns a reference to a variant setter object if the grid supports it.
 * Else returns a reference to an object that returns no supported variants
 * and throws an UnknownVariantError on set_variant.
 *
 * Attention: The returned reference points into "grid", thus, it is only
 * valid as long as "grid" is valid.
 */
VariantSetter &variants(grids::ParallelLCGrid *grid);

} // namespace repa