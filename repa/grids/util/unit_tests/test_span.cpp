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
#include <vector>

#include "../span.hpp"

TEST_CASE("span range")
{
    const std::vector<int> v{1, 2, 3, 4};

    const auto sp = repa::util::make_span(v);

    CHECK(std::equal(v.begin(), v.end(), sp.begin()));
    CHECK(std::equal(sp.begin(), sp.end(), v.begin()));
}

TEST_CASE("span reverse range")
{
    const std::vector<int> v{1, 2, 3, 4};

    const auto sp = repa::util::make_span(v);

    CHECK(std::equal(v.rbegin(), v.rend(), sp.rbegin()));
    CHECK(std::equal(sp.rbegin(), sp.rend(), v.rbegin()));
}

TEST_CASE("span data")
{
    const std::vector<int> v{1, 2, 3, 4};

    const auto sp = repa::util::make_span(v);

    CHECK(sp.data() == v.data());
    CHECK(sp.size() == v.size());

    for (size_t i = 0; i < v.size(); ++i) {
        CHECK(v[i] == sp[i]);
        CHECK_EQ(&v[i], &sp[i]);
    }
}

TEST_CASE("span empty")
{
    std::vector<int> v{1, 2, 3, 4};
    {
        const auto sp = repa::util::make_span(v);
        CHECK(sp.empty() == v.empty());
    }

    v.clear();
    {
        const auto sp = repa::util::make_span(v);
        CHECK(sp.empty() == v.empty());
        CHECK(sp.size() == v.size());
        CHECK(std::equal(v.begin(), v.end(), sp.begin()));
        CHECK(std::equal(sp.begin(), sp.end(), v.begin()));
        CHECK(std::equal(v.rbegin(), v.rend(), sp.rbegin()));
        CHECK(std::equal(sp.rbegin(), sp.rend(), v.rbegin()));
    }
}