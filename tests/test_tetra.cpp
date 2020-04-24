/**
 * Copyright 2017-2020 Steffen Hirschmann
 * Copyright 2020 Benjamin Vier
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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tetra
#include <boost/test/unit_test.hpp>

#include <array>
#include <chrono>
#include <random>

#include <repa/grids/util/tetra.hpp>

using namespace repa::util;
using repa::Vec3d;
using std::array;
typedef array<Vec3d, 8> octaVertices;

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
array<int, domains + 1>
ninsideDomains(array<octaVertices, domains> corners, int N, bool add)
{
    array<tetra::Octagon, domains> octs = {};
    for (int i = 0; i < domains; i++) {
        octs[i] = tetra::Octagon(corners[i]);
    }
    auto rnd = Randgen{};

    array<int, domains + 1> counter{0};
    for (int i = 0; i < N; ++i) {
        Vec3d p = {rnd(), rnd(), rnd()};
        if (add) {
            double factor = 1. / tetra::precision;
            p = {rnd() + factor, rnd() + factor, rnd() + factor};
        }
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
    PointArray();
    Vec3d randomPoint();
    octaVertices getVerticesAtPosition(int x, int y, int z);
};

Vec3d PointArray::randomPoint()
{
    // In the following, we assume that the midpoint (.5, .5, .5) of the domain
    // [0,1]^3 can be represented exactly in tetra internally.
    assert(std::floor(.5 * tetra::precision) == .5 * tetra::precision);

    Randgen rnd = Randgen{};
    Vec3d randVec = Vec3d{0, 0, 0};
    double sum = 0;
    for (int i = 0; i < 3; i++) {
        randVec[i] = 2 * rnd() - 1;
        sum += std::abs(randVec[i]);
    }
    double size = rnd();
    for (int i = 0; i < 3; i++) {
        // using l1-norm to prevent creating invalid Octagons.
        double norm1 = size * randVec[i] / sum;
        // cast from range (-1,1) to (0,1)
        randVec[i] = .5 + norm1 / 2;
        if (randVec[i] < .5) {
            randVec[i] += 1. / tetra::precision;
        }
    }
    return randVec;
}

PointArray::PointArray()
{
    for (int x = 0; x < size; x++) {
        for (int y = 0; y < size; y++) {
            for (int z = 0; z < size; z++) {
                point[x][y][z]
                    = Vec3d{{double(x) / 2, double(y) / 2, double(z) / 2}};
            }
        }
    }
    point[1][1][1] = randomPoint();
}

octaVertices PointArray::getVerticesAtPosition(int x, int y, int z)
{
    return {
        point[1 + x][1 + y][1 + z], point[0 + x][1 + y][1 + z],
        point[1 + x][0 + y][1 + z], point[0 + x][0 + y][1 + z],
        point[1 + x][1 + y][0 + z], point[0 + x][1 + y][0 + z],
        point[1 + x][0 + y][0 + z], point[0 + x][0 + y][0 + z],
    };
}

/**
 * Created one Octagon which covers 50% of a cube.
 * Test acceptance of N random points in this cube.
 * The number of accepted Points should be ~ 50%.
 */
BOOST_AUTO_TEST_CASE(test_tetra_1)
{

    tetra::Octagon r;
    BOOST_CHECK_EXCEPTION(r.contains({1., 2., 3.}), std::runtime_error,
                          empty_oct_message);

    // 50% of the volume of the unit cube
    octaVertices cs = {{{0., .5, 0.},
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
    array<octaVertices, 1> cArray{cs};
    auto acceptance = ninsideDomains<1>(cArray, N, false);
    double frac = static_cast<double>(acceptance[1]) / N;
    BOOST_CHECK((frac > .45 && frac < .55));
}

/**
 * Create two adjacent Octagons in a cube.
 * The cube is split in the first dimension.
 * The plane between the both Octagons is randomized.
 * Every random point in this cube should be accepted by exactly one Octagon.
 */
BOOST_AUTO_TEST_CASE(test_acceptance_of_two_domains_1)
{
    auto rnd = Randgen{};
    // p1-p4 define the randomized adjacent side of the Octagons
    Vec3d p1 = {rnd(), 1., 1.};
    Vec3d p2 = {rnd(), 0., 1.};
    Vec3d p3 = {rnd(), 1., 0.};
    Vec3d p4 = {rnd(), 0., 0.};
    // Create two Octagons
    array<octaVertices, 2> corners = {};
    corners[0] = {{{1., 1., 1.},
                   p1,
                   {1., 0., 1.},
                   p2,
                   {1., 1., 0.},
                   p3,
                   {1., 0., 0.},
                   p4}};
    corners[1] = {{p1,
                   {0., 1., 1.},
                   p2,
                   {0., 0., 1.},
                   p3,
                   {0., 1., 0.},
                   p4,
                   {0., 0., 0.}}};

    const int N = 10'000;
    array<int, 3> result = ninsideDomains<2>(corners, N, true);
    // All points should be accepted by exactly one Octagon.
    BOOST_CHECK(result[0] == 0);
    BOOST_CHECK(result[1] == N);
    BOOST_CHECK(result[2] == 0);
}

/**
 * Create two adjacent Octagons in a cube.
 * The cube is split in the second dimension.
 * The plane between the both Octagons is randomized.
 * Every random point in this cube should be accepted by exactly one Octagon.
 */
BOOST_AUTO_TEST_CASE(test_acceptance_of_two_domains_2)
{
    auto rnd = Randgen{};
    // p1-p4 define the randomized adjacent side of the Octagons
    Vec3d p1 = {1., rnd(), 1.};
    Vec3d p2 = {0., rnd(), 1.};
    Vec3d p3 = {1., rnd(), 0.};
    Vec3d p4 = {0., rnd(), 0.};
    // Create two Octagons
    array<octaVertices, 2> corners = {};
    corners[0] = {{{1., 1., 1.},
                   {0., 1., 1.},
                   p1,
                   p2,
                   {1., 1., 0.},
                   {0., 1., 0.},
                   p3,
                   p4}};
    corners[1] = {{p1,
                   p2,
                   {1., 0., 1.},
                   {0., 0., 1.},
                   p3,
                   p4,
                   {1., 0., 0.},
                   {0., 0., 0.}}};

    const int N = 10'000;
    array<int, 3> result = ninsideDomains<2>(corners, N, true);
    // All points should be accepted by exactly one Octagon.
    BOOST_CHECK(result[0] == 0);
    BOOST_CHECK(result[1] == N);
    BOOST_CHECK(result[2] == 0);
}

/**
 * Create two adjacent Octagons in a cube.
 * The cube is split in the third dimension.
 * The plane between the both Octagons is randomized.
 * Every random point in this cube should be accepted by exactly one Octagon.
 */
BOOST_AUTO_TEST_CASE(test_acceptance_of_two_domains_3)
{
    auto rnd = Randgen{};
    // p1-p4 define the randomized adjacent side of the Octagons
    Vec3d p1 = {1., 1., rnd()};
    Vec3d p2 = {0., 1., rnd()};
    Vec3d p3 = {1., 0., rnd()};
    Vec3d p4 = {0., 0., rnd()};
    // Create two Octagons
    array<octaVertices, 2> corners = {};
    corners[0] = {{{1., 1., 1.},
                   {0., 1., 1.},
                   {1., 0., 1.},
                   {0., 0., 1.},
                   p1,
                   p2,
                   p3,
                   p4}};
    corners[1] = {{p1,
                   p2,
                   p3,
                   p4,
                   {1., 1., 0.},
                   {0., 1., 0.},
                   {1., 0., 0.},
                   {0., 0., 0.}}};

    const int N = 10'000;
    array<int, 3> result = ninsideDomains<2>(corners, N, true);
    // All points should be accepted by exactly one Octagon.
    BOOST_CHECK(result[0] == 0);
    BOOST_CHECK(result[1] == N);
    BOOST_CHECK(result[2] == 0);
}

/**
 * Creating a cube filled with 8 Octagons.
 * Most of the points are randomized.
 * Every random point in this cube should be accepted by exactly one Octagon.
 */
BOOST_AUTO_TEST_CASE(test_tetra_3)
{
    using namespace repa;

    // Array with grid points of all Octagons.
    PointArray p{};

    // Create the 8 Octagons
    array<octaVertices, 8> corners = {};
    int index = 0;
    for (int x = 0; x < 2; x++) {
        for (int y = 0; y < 2; y++) {
            for (int z = 0; z < 2; z++) {
                corners[index] = p.getVerticesAtPosition(x, y, z);
                index++;
            }
        }
    }

    // Test acceptance of N points
    const int N = 10'000;
    array<int, 9> result = ninsideDomains<8>(corners, N, true);

    for (int i = 0; i < 9; i++) {
        // All points should be accepted exactly once.
        if (i == 1) {
            BOOST_CHECK(result[i] == N);
        }
        else {
            BOOST_CHECK(result[i] == 0);
        }
    }
}

/**
 * A Octagon should only accept points on 3 of its 6 sides.
 * These sides are predefined as the sides adjacent to the first vertex.
 */
BOOST_AUTO_TEST_CASE(test_tetra_4)
{
    octaVertices cs = {{{1., 1., 1.},
                        {0., 1., 1.},
                        {1., 0., 1.},
                        {0., 0., 1.},
                        {1., 1., 0.},
                        {0., 1., 0.},
                        {1., 0., 0.},
                        {0., 0., 0.}}};
    tetra::Octagon o = tetra::Octagon(cs);

    // These sides should NOT be accepted.
    BOOST_CHECK(!o.contains({0., .5, .5}));
    BOOST_CHECK(!o.contains({.5, 0., .5}));
    BOOST_CHECK(!o.contains({.5, .5, 0.}));
    // These sides should be accepted.
    BOOST_CHECK(o.contains({1., .5, .5}));
    BOOST_CHECK(o.contains({.5, 1., .5}));
    BOOST_CHECK(o.contains({.5, .5, 1.}));
}

/**
 * A Octagon should only accept points on 3 of its 6 sides.
 * These sides are predefined as the sides adjacent to the first vertex.
 */
BOOST_AUTO_TEST_CASE(test_validity_of_tetra)
{
    double max_cutoff = 2.;

    {
        // This Octagon should NOT be accepted.
        const octaVertices cs = {{{1., 1., 1.},
                                  {0., 1., 1.},
                                  {1., 0., 1.},
                                  {0., 0., 1.},
                                  {1., 1., 0.},
                                  {0., 1., 0.},
                                  {1., 0., 0.},
                                  {0., 0., 0.}}};
        BOOST_CHECK(!tetra::Octagon(cs, max_cutoff).is_valid());
    }
    {
        // This Octagon should NOT be accepted.
        const octaVertices cs2 = {{{12., 15., 15.},
                                   {0., 15., 15.},
                                   {15., 0., 15.},
                                   {0., 0., 15.},
                                   {15., 15., 0.},
                                   {0., 15., 0.},
                                   {15., 0., 0.},
                                   {0., 0., 0.}}};
        BOOST_CHECK(!tetra::Octagon(cs2, max_cutoff).is_valid());
    }
    {
        // This Octagon should be accepted.
        const octaVertices cs3 = {{{15., 15., 15.},
                                   {0., 15., 15.},
                                   {15., 0., 15.},
                                   {0., 0., 15.},
                                   {15., 15., 0.},
                                   {0., 15., 0.},
                                   {15., 0., 0.},
                                   {0., 0., 0.}}};
        BOOST_CHECK(tetra::Octagon(cs3, max_cutoff).is_valid());
    }
}

int main(int argc, char **argv)
{
    return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
