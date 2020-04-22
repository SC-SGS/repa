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

#include "grid_variants.hpp"

namespace {

/** Empty struct which allows not to set any variants.
 */
struct NoVariants : repa::VariantSetter {
    std::set<std::string> get_supported_variants() const
    {
        return {};
    };
    void set_variant(const std::string &)
    {
        throw repa::UnknownVariantError{};
    };
};

/** Meyer's singleton for NoVariants
 */
NoVariants &no_variants_singleton()
{
    static NoVariants novar;
    return novar;
}

} // namespace

namespace repa {

VariantSetter &variants(grids::ParallelLCGrid *grid)
{
    VariantSetter *p;
    if ((p = dynamic_cast<VariantSetter *>(grid)))
        return *p;
    else
        return no_variants_singleton();
}

} // namespace repa