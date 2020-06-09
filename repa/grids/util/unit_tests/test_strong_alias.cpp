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

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <doctest/doctest.h>
#include <sstream>

#include "../strong_alias.hpp"

struct _tag_type {
};

using s_alias = repa::util::StrongAlias<int, _tag_type>;

TEST_CASE("strong_alias uninitialized")
{
    const s_alias a;
#ifndef NDEBUG
    CHECK_FALSE(a.is_initialized());
#endif
}

TEST_CASE("strong_alias constuctors")
{
    const s_alias a{1};
#ifndef NDEBUG
    CHECK(a.is_initialized());
#endif
    CHECK(a == 1);
    CHECK(static_cast<int>(a) == 1);
    CHECK(a.value() == 1);

    const s_alias b{a};
#ifndef NDEBUG
    CHECK(b.is_initialized());
#endif
    CHECK(b == 1);
    CHECK(static_cast<int>(b) == 1);
    CHECK(b.value() == 1);
}

TEST_CASE("strong_alias value constructor")
{
    const s_alias a{UINT64_C(0)};
#ifndef NDEBUG
    CHECK(a.is_initialized());
#endif
    CHECK(a == 0);
    CHECK(static_cast<int>(a) == 0);
    CHECK(a.value() == 0);
}

TEST_CASE("strong_alias assignment")
{
    const s_alias a{1};
    s_alias b;

    b = a;
#ifndef NDEBUG
    CHECK(b.is_initialized());
#endif
    CHECK(b == 1);
    CHECK(static_cast<int>(a) == 1);
    CHECK(a.value() == 1);
}

TEST_CASE("strong_alias increment")
{
    s_alias a{1};

    CHECK_EQ(a++, 1);
    CHECK(a == 2);
    CHECK(static_cast<int>(a) == 2);
    CHECK(a.value() == 2);

    CHECK_EQ(++a, 3);
    CHECK(a == 3);
    CHECK(static_cast<int>(a) == 3);
    CHECK(a.value() == 3);
}

TEST_CASE("strong_alias assignment arith")
{
    s_alias a{1};

    CHECK_EQ(a += s_alias{2}, 3);
    CHECK(a == 3);
    CHECK(static_cast<int>(a) == 3);
    CHECK(a.value() == 3);

    CHECK_EQ(a -= s_alias{2}, 1);
    CHECK(a == 1);
    CHECK(static_cast<int>(a) == 1);
    CHECK(a.value() == 1);
}

TEST_CASE("strong_alias cmp")
{
    const s_alias a{1};

    const std::vector<int> ints{-1, 0, 1, 2, 3};

    for (const int i : ints) {
        CHECK_EQ(a == s_alias{i}, a == i);
        CHECK_EQ(a != s_alias{i}, a != i);
        CHECK_EQ(a < s_alias{i}, a < i);
        CHECK_EQ(a <= s_alias{i}, a <= i);
        CHECK_EQ(a > s_alias{i}, a > i);
        CHECK_EQ(a >= s_alias{i}, a >= i);
    }
}

TEST_CASE("strong_alias serialization")
{
    const s_alias a{12};

    // Serialization
    std::stringstream sstr;
    {
        boost::archive::text_oarchive oar{sstr};
        oar << a;
    }

    s_alias b;
    {
        boost::archive::text_iarchive iar{sstr};
        iar >> b;
    }

#ifndef NDEBUG
    CHECK(b.is_initialized());
#endif
    CHECK(a == b);
}