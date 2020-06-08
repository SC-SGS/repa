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

#include "../ioptional.hpp"

using optional_int = repa::util::ioptional<repa::rank_type>;

TEST_CASE("ioptional empty")
{
    auto io = optional_int{};
    CHECK_FALSE(io.has_value());
    io = 12;
    CHECK(io.has_value());

    auto io2 = optional_int{12};
    CHECK(io2.has_value());
}

TEST_CASE("ioptional value")
{
    auto io = optional_int{};
    CHECK(io.value_or(-255) == -255);

    CHECK_THROWS_WITH_AS((io.value_or_throw<std::runtime_error>("test")),
                         "test", std::runtime_error);

    io = 12;
    CHECK(io.value() == 12);
    CHECK(io.value_or(-255) == 12);
    CHECK(io.ensure_value() == 12);
    CHECK(io.value_or_throw<std::runtime_error>("test") == 12);
}

TEST_CASE("ioptional operator bool")
{
    auto io = optional_int{};
    CHECK_FALSE(io);
    io = 12;
    CHECK(io);
}

TEST_CASE("ioptional deref")
{
    auto io = optional_int{12};
    CHECK_EQ(*io, 12);
}

TEST_CASE("ioptional reset")
{
    auto io = optional_int{12};
    CHECK(io.has_value());
    io.reset();
    CHECK_FALSE(io.has_value());
}

TEST_CASE("ioptional cmp")
{
    auto io = optional_int{};
    auto io2 = optional_int{12};
    CHECK(io != io2);

    io = io2;
    CHECK(io == io2);
}