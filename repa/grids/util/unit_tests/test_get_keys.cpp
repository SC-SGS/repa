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
#include <map>
#include <string>
#include <unordered_map>

#include "../get_keys.hpp"

TEST_CASE("get_keys map")
{
    const std::map<int, std::string> m{{1, "one"}, {2, "two"}, {3, "three"}};
    auto keys = repa::util::get_keys(m);
    for (const auto &k : m) {
        CHECK(keys.find(k.first) != std::end(keys));
        keys.erase(k.first);
    }
    CHECK(keys.empty());
}

TEST_CASE("get_keys unordered_map")
{
    const std::unordered_map<int, std::string> m{
        {1, "one"}, {2, "two"}, {3, "three"}};
    auto keys = repa::util::get_keys(m);
    for (const auto &k : m) {
        CHECK(keys.find(k.first) != std::end(keys));
        keys.erase(k.first);
    }
    CHECK(keys.empty());
}
