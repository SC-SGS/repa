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

#include <doctest/doctest.h>

#include "../linearize.hpp"

TEST_CASE("linearize bijectivity")
{
    const auto grid = repa::Vec3i{10, 10, 10};
    std::vector<bool> hit(1000, false);

    repa::Vec3i x;
    for (x[0] = 0; x[0] < grid[0]; x[0]++) {
        for (x[1] = 0; x[1] < grid[1]; x[1]++) {
            for (x[2] = 0; x[2] < grid[2]; x[2]++) {
                auto lin_i = repa::util::linearize(x, grid);
                CHECK(lin_i >= 0);
                CHECK(lin_i < 1000);
                CHECK_FALSE(hit[lin_i]);
                hit[lin_i] = true;
            }
        }
    }
}

TEST_CASE("linearize/unlinearize inverse")
{
    const auto grid = repa::Vec3i{10, 10, 10};
    std::vector<bool> used(1000, false);

    repa::Vec3i x;
    for (x[0] = 0; x[0] < grid[0]; x[0]++) {
        for (x[1] = 0; x[1] < grid[1]; x[1]++) {
            for (x[2] = 0; x[2] < grid[2]; x[2]++) {
                auto y = repa::util::unlinearize(repa::util::linearize(x, grid),
                                                 grid);
                CHECK(x == y);
            }
        }
    }
}