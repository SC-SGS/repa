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
#include <set>

#include "../set_union.hpp"

TEST_CASE("set_union superset_eq")
{
    std::set<int> s1{1, 2, 3, 4, 5}, s2{4, 5, 6, 7, 8, 9};
    const auto su = repa::util::set_union(s1, s2);

    for (const auto i : s1)
        CHECK(su.find(i) != su.end());
    for (const auto i : s2)
        CHECK(su.find(i) != su.end());
}

TEST_CASE("set_union subset_eq")
{
    std::set<int> s1{1, 2, 3, 4, 5}, s2{4, 5, 6, 7, 8, 9};
    auto su = repa::util::set_union(s1, s2);

    for (const auto i : s1)
        su.erase(i);
    for (const auto i : s2)
        su.erase(i);
    CHECK(su.empty());
}