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
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE strong_alias
#include <boost/test/unit_test.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <sstream>

#include <repa/grids/util/strong_alias.hpp>

struct _tag_type {
};

using s_alias = repa::util::StrongAlias<int, _tag_type>;

BOOST_AUTO_TEST_CASE(uninitialized)
{
#ifndef NDEBUG
    const s_alias a;
    BOOST_TEST((!a.is_initialized()));
#endif
}

BOOST_AUTO_TEST_CASE(constuctors)
{
    const s_alias a{1};
#ifndef NDEBUG
    BOOST_TEST(a.is_initialized());
#endif
    BOOST_TEST(a == 1);
    BOOST_TEST(static_cast<int>(a) == 1);
    BOOST_TEST(a.value() == 1);

    const s_alias b{a};
#ifndef NDEBUG
    BOOST_TEST(b.is_initialized());
#endif
    BOOST_TEST(b == 1);
    BOOST_TEST(static_cast<int>(b) == 1);
    BOOST_TEST(b.value() == 1);
}

BOOST_AUTO_TEST_CASE(value_constructor)
{
    const s_alias a{UINT64_C(0)};
#ifndef NDEBUG
    BOOST_TEST(a.is_initialized());
#endif
    BOOST_TEST(a == 0);
    BOOST_TEST(static_cast<int>(a) == 0);
    BOOST_TEST(a.value() == 0);
}

BOOST_AUTO_TEST_CASE(assignment)
{
    const s_alias a{1};
    s_alias b;

    b = a;
#ifndef NDEBUG
    BOOST_TEST(b.is_initialized());
#endif
    BOOST_TEST(b == 1);
    BOOST_TEST(static_cast<int>(a) == 1);
    BOOST_TEST(a.value() == 1);
}

BOOST_AUTO_TEST_CASE(increment)
{
    s_alias a{1};

    BOOST_CHECK_EQUAL(a++, 1);
    BOOST_TEST(a == 2);
    BOOST_TEST(static_cast<int>(a) == 2);
    BOOST_TEST(a.value() == 2);

    BOOST_CHECK_EQUAL(++a, 3);
    BOOST_TEST(a == 3);
    BOOST_TEST(static_cast<int>(a) == 3);
    BOOST_TEST(a.value() == 3);
}

BOOST_AUTO_TEST_CASE(assignment_arith)
{
    s_alias a{1};

    BOOST_CHECK_EQUAL(a += s_alias{2}, 3);
    BOOST_TEST(a == 3);
    BOOST_TEST(static_cast<int>(a) == 3);
    BOOST_TEST(a.value() == 3);

    BOOST_CHECK_EQUAL(a -= s_alias{2}, 1);
    BOOST_TEST(a == 1);
    BOOST_TEST(static_cast<int>(a) == 1);
    BOOST_TEST(a.value() == 1);
}

BOOST_AUTO_TEST_CASE(cmp)
{
    const s_alias a{1};

    const std::vector<int> ints{-1, 0, 1, 2, 3};

    for (const int i : ints) {
        BOOST_CHECK_EQUAL(a == s_alias{i}, a == i);
        BOOST_CHECK_EQUAL(a != s_alias{i}, a != i);
        BOOST_CHECK_EQUAL(a < s_alias{i}, a < i);
        BOOST_CHECK_EQUAL(a <= s_alias{i}, a <= i);
        BOOST_CHECK_EQUAL(a > s_alias{i}, a > i);
        BOOST_CHECK_EQUAL(a >= s_alias{i}, a >= i);
    }
}

BOOST_AUTO_TEST_CASE(serialization)
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
    BOOST_TEST(b.is_initialized());
#endif
    BOOST_TEST(a == b);
}