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
using repa::Vec3d;
using std::array;

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

template <int domains>
array<int, domains + 1> ninsideDomains(array<Vec3d, 8> corners[domains], int N)
{
    array<tetra::Octagon, domains> octs = {};
    for (int i = 0; i < domains; i++) {
        octs[i] = tetra::Octagon(corners[i]);
    }
    auto rnd = Randgen{};

    array<int, domains + 1> counter{0};
    for (int i = 0; i < N; ++i) {
        Vec3d p = {rnd(), rnd(), rnd()};
        int count = 0;
        for (int k = 0; k < domains; ++k) {
            if (octs[k].contains(p)) {
                count++;
            }
        }
        counter[count]++;
    }
    return counter;
}

struct PointArray {
    static const int size = 3;
    array<array<array<Vec3d, size>, size>, size> point = {{{}}};
    Randgen rnd = Randgen{};
    PointArray();
    double valueFor(int i)
    {
        if (i == 1) {
            return 0.25 + 0.5 * rnd();
        }
        return double(i) / 2; // map i={0|2} to return {0|1}
    };
    array<Vec3d, 8> getVerticesAtPosition(int x, int y, int z);
};

PointArray::PointArray()
{
    for (int x = 0; x < size; x++) {
        for (int y = 0; y < size; y++) {
            for (int z = 0; z < size; z++) {
                point[x][y][z] = Vec3d{{valueFor(x), valueFor(y), valueFor(z)}};
            }
        }
    }
    point[1][1][1] = Vec3d{{0.5, 0.5, 0.5}};
}

array<Vec3d, 8> PointArray::getVerticesAtPosition(int x, int y, int z)
{
    return {
        point[0 + x][0 + y][0 + z], point[1 + x][0 + y][0 + z],
        point[0 + x][1 + y][0 + z], point[1 + x][1 + y][0 + z],
        point[0 + x][0 + y][1 + z], point[1 + x][0 + y][1 + z],
        point[0 + x][1 + y][1 + z], point[1 + x][1 + y][1 + z],
    };
};

BOOST_AUTO_TEST_CASE(test_tetra_1)
{

    tetra::Octagon r;
    BOOST_CHECK_EXCEPTION(r.contains({1., 2., 3.}), std::runtime_error,
                          empty_oct_message);

    // 50% of the volume of the unit cube
    array<Vec3d, 8> cs = {{{0., .5, 0.},
                           {0., 1., .5},
                           {0., 0., .5},
                           {0., .5, 1.},
                           {1., .5, 0.},
                           {1., 1., .5},
                           {1., 0., .5},
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
    array<Vec3d, 8> cArray[1]{cs};
    auto acceptance = ninsideDomains<1>(cArray, N);
    double frac = static_cast<double>(acceptance[1]) / N;
    BOOST_CHECK((frac > .45 && frac < .55));
}

BOOST_AUTO_TEST_CASE(test_tetra_2)
{
    auto rnd = Randgen{};
    Vec3d p1 = {rnd(), 0., 0.};
    Vec3d p2 = {rnd(), 1., 0.};
    Vec3d p3 = {rnd(), 0., 1.};
    Vec3d p4 = {rnd(), 1., 1.};
    array<Vec3d, 8> corners[2] = {};
    corners[0] = {{{0., 0., 0.},
                   p1,
                   {0., 1., 0.},
                   p2,
                   {0., 0., 1.},
                   p3,
                   {0., 1., 1.},
                   p4}};
    corners[1] = {{p1,
                   {1., 0., 0.},
                   p2,
                   {1., 1., 0.},
                   p3,
                   {1., 0., 1.},
                   p4,
                   {1., 1., 1.}}};

    const int N = 10'000;
    array<int, 3> result = ninsideDomains<2>(corners, N);
    BOOST_CHECK(result[0] == 0);
    BOOST_CHECK(result[1] == N);
    BOOST_CHECK(result[2] == 0);
}

BOOST_AUTO_TEST_CASE(test_tetra_3)
{
    using namespace repa;
    auto rnd = Randgen{};

    PointArray p{};

    array<Vec3d, 8> corners[8] = {};
    int index = 0;
    for (int x = 0; x < 2; x++) {
        for (int y = 0; y < 2; y++) {
            for (int z = 0; z < 2; z++) {
                corners[index] = p.getVerticesAtPosition(x, y, z);
                index++;
            }
        }
    }

    const int N = 10'000;
    array<int, 9> result = ninsideDomains<8>(corners, N);
    BOOST_CHECK(result[0] == 0);
    BOOST_CHECK(result[1] == N);

    for (int i = 2; i < 9; i++) {
        BOOST_CHECK(result[i] == 0);
    }
}

BOOST_AUTO_TEST_CASE(test_tetra_4)
{
    array<Vec3d, 8> cs = {{{0., 0., 0.},
                           {1., 0., 0.},
                           {0., 1., 0.},
                           {1., 1., 0.},
                           {0., 0., 1.},
                           {1., 0., 1.},
                           {0., 1., 1.},
                           {1., 1., 1.}}};
    tetra::Octagon o = tetra::Octagon(cs);
    BOOST_CHECK(o.contains({0., .5, .5}));
    BOOST_CHECK(o.contains({.5, 0., .5}));
    BOOST_CHECK(o.contains({.5, .5, 0.}));
    BOOST_CHECK(!o.contains({1., .5, .5}));
    BOOST_CHECK(!o.contains({.5, 1., .5}));
    BOOST_CHECK(!o.contains({.5, .5, 1.}));
}