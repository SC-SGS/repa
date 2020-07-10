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
#define BOOST_TEST_MODULE span
#include <boost/test/unit_test.hpp>

#include <vector>

#include <repa/grids/util/span.hpp>

BOOST_AUTO_TEST_CASE(range)
{
    const std::vector<int> v{1, 2, 3, 4};

    const auto sp = repa::util::make_span(v);

    BOOST_TEST(std::equal(v.begin(), v.end(), sp.begin()));
    BOOST_TEST(std::equal(sp.begin(), sp.end(), v.begin()));
}

BOOST_AUTO_TEST_CASE(reverse_range)
{
    const std::vector<int> v{1, 2, 3, 4};

    const auto sp = repa::util::make_span(v);

    BOOST_TEST(std::equal(v.rbegin(), v.rend(), sp.rbegin()));
    BOOST_TEST(std::equal(sp.rbegin(), sp.rend(), v.rbegin()));
}

BOOST_AUTO_TEST_CASE(data)
{
    const std::vector<int> v{1, 2, 3, 4};

    const auto sp = repa::util::make_span(v);

    BOOST_TEST(sp.data() == v.data());
    BOOST_TEST(sp.size() == v.size());

    for (size_t i = 0; i < v.size(); ++i) {
        BOOST_TEST(v[i] == sp[i]);
        BOOST_CHECK_EQUAL(&v[i], &sp[i]);
    }
}

BOOST_AUTO_TEST_CASE(empty)
{
    std::vector<int> v{1, 2, 3, 4};
    {
        const auto sp = repa::util::make_span(v);
        BOOST_TEST(sp.empty() == v.empty());
    }

    v.clear();
    {
        const auto sp = repa::util::make_span(v);
        BOOST_TEST(sp.empty() == v.empty());
        BOOST_TEST(sp.size() == v.size());
        BOOST_TEST(std::equal(v.begin(), v.end(), sp.begin()));
        BOOST_TEST(std::equal(sp.begin(), sp.end(), v.begin()));
        BOOST_TEST(std::equal(v.rbegin(), v.rend(), sp.rbegin()));
        BOOST_TEST(std::equal(sp.rbegin(), sp.rend(), v.rbegin()));
    }
}