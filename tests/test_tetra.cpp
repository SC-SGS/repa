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
 * Checks tetra.[ch]pp
 */

#define BOOST_TEST_MODULE tetra

#include <array>
#include <chrono>
#include <random>

#include <boost/test/included/unit_test.hpp>
#include <repa/grids/util/tetra.hpp>

using namespace repa::util;

bool empty_oct_message(const std::runtime_error &e)
{
    return std::string{e.what()} == "contains() on empty octagon";
}

struct Randgen {
    Randgen()
        : mt(std::chrono::high_resolution_clock::now()
                 .time_since_epoch()
                 .count()),
          d(0, 1)
    {
    }
    double operator()()
    {
        return d(mt);
    }

private:
    std::mt19937 mt;
    std::uniform_real_distribution<> d;
};

int ninside(const tetra::Octagon &o, int N)
{
    auto rnd = Randgen{};
    int n = 0;

    for (int i = 0; i < N; ++i) {
        if (o.contains({rnd(), rnd(), rnd()}))
            n++;
    }
    return n;
}

BOOST_AUTO_TEST_CASE(test_tetra)
{

    tetra::Octagon r;
    BOOST_CHECK_EXCEPTION(r.contains({1., 2., 3.}), std::runtime_error,
                          empty_oct_message);

    // 50% of the volume of the unit cube
    std::array<repa::Vec3d, 8> cs = {{{0., .5, 0.},
                                      {0., 0., .5},
                                      {0., 1., .5},
                                      {0., .5, 1.},
                                      {1., .5, 0.},
                                      {1., 0., .5},
                                      {1., 1., .5},
                                      {1., .5, 1.}}};
    auto o = tetra::Octagon{cs};

    BOOST_CHECK(o.contains({.5, .5, .5}));

    BOOST_CHECK(!o.contains({.2, .2, .2}));
    BOOST_CHECK(!o.contains({.2, .2, .8}));
    BOOST_CHECK(!o.contains({.2, .8, .2}));
    BOOST_CHECK(!o.contains({.2, .8, .8}));
    BOOST_CHECK(!o.contains({.8, .2, .2}));
    BOOST_CHECK(!o.contains({.8, .2, .8}));
    BOOST_CHECK(!o.contains({.8, .8, .2}));
    BOOST_CHECK(!o.contains({.8, .2, .8}));
    BOOST_CHECK(!o.contains({.8, .8, .2}));
    BOOST_CHECK(!o.contains({.8, .8, .8}));

    const int N = 10'000;
    double frac = static_cast<double>(ninside(o, N)) / N;
    BOOST_CHECK((frac > .45 && frac < .55));
}
