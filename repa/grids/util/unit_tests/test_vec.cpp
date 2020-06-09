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

#include <array>
#include <chrono>
#include <random>
#include <sstream>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "common_types.hpp"

TEST_CASE("vec size")
{
    repa::Vec3d v;

    CHECK(v.size() == 3);
    CHECK(v.max_size() == 3);
    CHECK(!v.empty());
    CHECK(v.data() != nullptr);

    CHECK(static_cast<size_t>(std::distance(v.begin(), v.end())) == v.size());
    CHECK(std::distance(v.begin(), v.end())
          == static_cast<ptrdiff_t>(v.size()));
}

TEST_CASE("vec iterators")
{
    repa::Vec3i v{1, 2, 3};

    int i = 1;
    for (const auto x : v) {
        CHECK(x == i);
        i++;
    }

    for (repa::Vec3i::reverse_iterator it = v.rbegin(); it != v.rend(); it++) {
        i--;
        CHECK(*it == i);
    }

    i = 1;
    const repa::Vec3i &vconst = v;
    for (const auto x : vconst) {
        CHECK(x == i);
        i++;
    }

    for (repa::Vec3i::const_reverse_iterator it = vconst.rbegin();
         it != vconst.rend(); it++) {
        i--;
        CHECK(*it == i);
    }
}

TEST_CASE("vec uninitialized")
{
    repa::Vec3d v;

    for (const auto x : v)
        CHECK(x == 0.0);

    const repa::Vec3d &vconst = v;
    for (const auto x : vconst)
        CHECK(x == 0.0);
}

TEST_CASE("vec memory loc")
{
    repa::Vec3d v;
    repa::Vec3d w{v};
    w[1] = 1.0;

    const repa::Vec3d &wconst = w;

    for (size_t i = 0; i < w.size(); ++i) {
        CHECK_EQ(&w.data()[i], &w[i]);
        CHECK_NE(&w[i], &v[i]);
        CHECK_EQ(&wconst.data()[i], &w[i]);
        CHECK_EQ(&w.as_array().data()[i], &w[i]);
        CHECK_EQ(&wconst.as_array().data()[i], &w[i]);
    }
}

TEST_CASE("vec cmp")
{
    repa::Vec3d v;
    repa::Vec3d w{v};

    CHECK(v == v);
    CHECK(w == w);
    CHECK(v == w);

    for (size_t i = 0; i < w.size(); ++i) {
        w[i] = 1.0;
        CHECK(v != w);

        w[i] = v[i];
        CHECK(v == w);
    }
}

TEST_CASE("vec move")
{
    repa::Vec3d w;
    w[1] = 1.0;

    repa::Vec3d w_copy{w};

    // Move constructor
    repa::Vec3d u{std::move(w)};
    CHECK(u[0] == 0.0);
    CHECK(u[1] == 1.0);
    CHECK(u[2] == 0.0);

    CHECK(u == w_copy);
}

TEST_CASE("vec assignment")
{
    repa::Vec3i i{1, 2, 3};
    // Assignment operator
    repa::Vec3i k;
    k = i;
    CHECK(i == k);
    k[0] -= 1;
    CHECK(i != k);
}

TEST_CASE("vec serialization")
{
    repa::Vec3i k{1, 2, 3};
    std::stringstream sstr;
    {
        boost::archive::text_oarchive oar{sstr};
        oar << k;
    }

    repa::Vec3i k2;
    CHECK(k != k2);
    {
        boost::archive::text_iarchive iar{sstr};
        iar >> k2;
    }

    CHECK(k == k2);
}
