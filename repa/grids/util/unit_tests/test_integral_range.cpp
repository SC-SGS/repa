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

#include "common_types.hpp"
#include <doctest/doctest.h>

TEST_CASE("integral_range domain")
{
    for (repa::fs_neighidx::value_type i = 0; i <= 26; ++i) {
        repa::fs_neighidx ni{i};
        CHECK(static_cast<repa::fs_neighidx::value_type>(ni) == i);
    }

#ifndef NDEBUG
    repa::fs_neighidx ni = 0;
    CHECK_THROWS_AS(ni = 27, std::domain_error);
    CHECK_THROWS_AS(ni = -1, std::domain_error);
#endif
}
