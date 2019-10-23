/**
 * Copyright 2017-2019 Steffen Hirschmann
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

/**
 * Checks Vec type in common_types.hpp
 */

#define BOOST_TEST_MODULE vec

#include <array>
#include <chrono>
#include <random>
#include <sstream>

#include <boost/test/included/unit_test.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <repa/repa.hpp>

BOOST_AUTO_TEST_CASE(test_vec)
{
    repa::Vec3d v;

    // Size
    BOOST_CHECK(v.size() == 3);
    BOOST_CHECK(v.max_size() == 3);
    BOOST_CHECK(!v.empty());
    BOOST_CHECK(v.data() != nullptr);

    // Iterators
    for (const auto x : v)
        BOOST_CHECK(x == 0.0);

    const repa::Vec3d &vconst = v;
    for (const auto x : vconst)
        BOOST_CHECK(x == 0.0);

    repa::Vec3d w{v};
    w[1] = 1.0;

    // Reverse iterators
    for (repa::Vec3d::reverse_iterator i = v.rbegin(); i != v.rend(); i++)
        BOOST_CHECK(*i == 0.0);
    for (repa::Vec3d::const_reverse_iterator i = v.crbegin(); i != v.crend();
         i++)
        BOOST_CHECK(*i == 0.0);

    // Various methods to get the memory locations
    const repa::Vec3d &wconst = w;
    for (size_t i = 0; i < w.size(); ++i) {
        BOOST_CHECK((&w.data()[i] == &w[i]) && (&w[i] != &v[i]));
        BOOST_CHECK(&wconst.data()[i] == &w[i]);
        BOOST_CHECK(&w.as_array().data()[i] == &w[i]);
        BOOST_CHECK(&wconst.as_array().data()[i] == &w[i]);
    }

    BOOST_CHECK(std::distance(w.begin(), w.end()) == w.size());

    // Comparisons
    BOOST_CHECK(v == v);
    BOOST_CHECK(w == w);
    BOOST_CHECK(v != w);
    BOOST_CHECK((w > v) && (w >= v));
    BOOST_CHECK(!(w < v) && !(w <= v));

    // Move constructor
    repa::Vec3d u{std::move(w)};
    BOOST_CHECK(u[0] == 0.0);
    BOOST_CHECK(u[1] == 1.0);
    BOOST_CHECK(u[2] == 0.0);

    // Int Vecs and copy constructor
    repa::Vec3i i{1, 2, 3};
    BOOST_CHECK(i[0] == 1);
    BOOST_CHECK(i[1] == 2);
    BOOST_CHECK(i[2] == 3);
    repa::Vec3i j{i};
    BOOST_CHECK((i == j) && !(i != j));
    j[0] += 1;
    BOOST_CHECK((i != j) && !(i == j));

    // Assignment operator
    repa::Vec3i k;
    k = i;
    BOOST_CHECK((i == k));
    k[0] -= 1;
    BOOST_CHECK((i != k));

    // Serialization
    std::stringstream sstr;
    {
        boost::archive::text_oarchive oar{sstr};
        oar << k;
    }

    repa::Vec3i k2;
    BOOST_CHECK((k2[0] == 0) && (k2[1] == 0) && (k2[2] == 0));
    BOOST_CHECK(k != k2);
    {
        boost::archive::text_iarchive iar{sstr};
        iar >> k2;
    }

    BOOST_CHECK((k == k2));
}
