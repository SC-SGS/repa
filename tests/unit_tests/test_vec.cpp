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
#define BOOST_TEST_MODULE vec
#include <boost/test/unit_test.hpp>

#include <array>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <chrono>
#include <random>
#include <sstream>

#include <repa/common_types.hpp>

BOOST_AUTO_TEST_CASE(size)
{
    repa::Vec3d v;

    BOOST_TEST(v.size() == 3);
    BOOST_TEST(v.max_size() == 3);
    BOOST_TEST((!v.empty()));
    BOOST_TEST(v.data() != nullptr);

    BOOST_TEST(static_cast<size_t>(std::distance(v.begin(), v.end()))
               == v.size());
    BOOST_TEST(std::distance(v.begin(), v.end())
               == static_cast<ptrdiff_t>(v.size()));
}

BOOST_AUTO_TEST_CASE(iterators)
{
    repa::Vec3i v{1, 2, 3};

    int i = 1;
    for (const auto x : v) {
        BOOST_TEST(x == i);
        i++;
    }

    for (repa::Vec3i::reverse_iterator it = v.rbegin(); it != v.rend(); it++) {
        i--;
        BOOST_TEST(*it == i);
    }

    i = 1;
    const repa::Vec3i &vconst = v;
    for (const auto x : vconst) {
        BOOST_TEST(x == i);
        i++;
    }

    for (repa::Vec3i::const_reverse_iterator it = vconst.rbegin();
         it != vconst.rend(); it++) {
        i--;
        BOOST_TEST(*it == i);
    }
}

BOOST_AUTO_TEST_CASE(uninitialized)
{
    repa::Vec3d v;

    for (const auto x : v)
        BOOST_TEST(x == 0.0);

    const repa::Vec3d &vconst = v;
    for (const auto x : vconst)
        BOOST_TEST(x == 0.0);
}

BOOST_AUTO_TEST_CASE(memory_loc)
{
    repa::Vec3d v;
    repa::Vec3d w{v};
    w[1] = 1.0;

    const repa::Vec3d &wconst = w;

    for (size_t i = 0; i < w.size(); ++i) {
        BOOST_CHECK_EQUAL(&w.data()[i], &w[i]);
        BOOST_CHECK_NE(&w[i], &v[i]);
        BOOST_CHECK_EQUAL(&wconst.data()[i], &w[i]);
        BOOST_CHECK_EQUAL(&w.as_array().data()[i], &w[i]);
        BOOST_CHECK_EQUAL(&wconst.as_array().data()[i], &w[i]);
    }
}

BOOST_AUTO_TEST_CASE(cmp)
{
    repa::Vec3d v;
    repa::Vec3d w{v};

    BOOST_TEST(v == v);
    BOOST_TEST(w == w);
    BOOST_TEST(v == w);

    for (size_t i = 0; i < w.size(); ++i) {
        w[i] = 1.0;
        BOOST_TEST(v != w);

        w[i] = v[i];
        BOOST_TEST(v == w);
    }
}

BOOST_AUTO_TEST_CASE(move)
{
    repa::Vec3d w;
    w[1] = 1.0;

    repa::Vec3d w_copy{w};

    // Move constructor
    repa::Vec3d u{std::move(w)};
    BOOST_TEST(u[0] == 0.0);
    BOOST_TEST(u[1] == 1.0);
    BOOST_TEST(u[2] == 0.0);

    BOOST_TEST(u == w_copy);
}

BOOST_AUTO_TEST_CASE(assignment)
{
    repa::Vec3i i{1, 2, 3};
    // Assignment operator
    repa::Vec3i k;
    k = i;
    BOOST_TEST(i == k);
    k[0] -= 1;
    BOOST_TEST(i != k);
}

BOOST_AUTO_TEST_CASE(serialization)
{
    repa::Vec3i k{1, 2, 3};
    std::stringstream sstr;
    {
        boost::archive::text_oarchive oar{sstr};
        oar << k;
    }

    repa::Vec3i k2;
    BOOST_TEST(k != k2);
    {
        boost::archive::text_iarchive iar{sstr};
        iar >> k2;
    }

    BOOST_TEST(k == k2);
}
