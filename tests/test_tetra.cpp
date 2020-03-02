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

std::array<int, 3> ninside2Domains(std::array<, 8> corner1, std::array<, 8> corner2, int N)
{
    auto o1 = tetra::Octagon(corner1);
    auto o2 = tetra::Octagon(corner2);
    auto rnd = Randgen{};

    std::array<int, 3> counter{ {0,0,0} };

    for (int i = 0; i < N; ++i) {
         p = { rnd(), rnd(), rnd() };
        bool c1 = o1.contains(p);
        bool c2 = o2.contains(p);
        counter[c1 + c2]++;
    }
    return counter;
}

std::array<int, 9> ninside8Domains(std::array<, 8> corners[8], int N)
{
    std::array<tetra::Octagon, 8> octs = {};
    for (int i = 0; i < 8; i++) {
        octs[i] = tetra::Octagon(corners[i]);
    }
    auto rnd = Randgen{};

    std::array<int, 9> counter{ {0,0,0,0,0,0,0,0,0} };
    std::ofstream csv("doublePoints.csv");
    for (int i = 0; i < N; ++i) {
         p = { rnd(), rnd(), rnd() };
        int count = 0;
        for (int k = 0; k < 8; ++k) {
            if (octs[k].contains(p)) { count++; }
        }
        counter[count]++;
        if (count != 1) {
            csv << std::to_string(p[0]) << "," << std::to_string(p[1]) << "," << std::to_string(p[2]) << "\n";
        }
    }
    csv.close();
    return counter;
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

BOOST_AUTO_TEST_CASE(test_tetra)
{
    auto rnd = Randgen{};
    auto p1 = { rnd(), 0, 0 };
    auto p2 = { rnd(), 1, 0 };
    auto p3 = { rnd(), 0, 1 };
    auto p4 = { rnd(), 1, 1 };
    std::array<repa::Vec3d, 8> corner1 = { {0, 0, 0},
                                p1,
                                {0, 1, 0},
                                p2,
                                {0, 0, 1},
                                p3,
                                {0, 1, 1},
                                p4};
    std::array<repa::Vec3d, 8> corner2 = { p1,
                                {1, 0, 0},
                                p2,
                                {1, 1, 0},
                                p3,
                                {1, 0, 1},
                                p4,
                                {1, 1, 1}};

    const int N = 10'000;
    std::array<int, 3> result = ninsideTwoDomains(corner1, corner2, N);
    BOOST_CHECK(result[0] == 0);
    BOOST_CHECK(result[1] == N);
    BOOST_CHECK(result[2] == 0);
}

BOOST_AUTO_TEST_CASE(test_tetra)
{
    using namespace repa;
    auto rnd = Randgen{};

    Vec3d point[3][3][3]{};
    for (int i = 0; i < 3; i++) {
        for (int k = 0; k < 3; k++) {
            for (int j = 0; j < 3; j++) {
                point[i][k][j] = Vec3d{ double(i) / 2.0,double(k) / 2.0,double(j) / 2.0 };
            }
        }
    }
    point[1][1][1] = Vec3d{ rnd(), rnd(), rnd() };


    std::array<, 8> corners[8] = {};
    corners[0] = { point[0][0][0], point[1][0][0], point[0][1][0], point[1][1][0],
        point[0][0][1], point[1][0][1], point[0][1][1], point[1][1][1] };
    corners[1] = { point[1][0][0], point[2][0][0], point[1][1][0], point[2][1][0],
        point[1][0][1], point[2][0][1], point[1][1][1], point[2][1][1] };
    corners[2] = { point[0][1][0], point[1][1][0], point[0][2][0], point[1][2][0],
        point[0][1][1], point[1][1][1], point[0][2][1], point[1][2][1] };
    corners[3] = { point[1][1][0], point[2][1][0], point[1][2][0], point[2][2][0],
        point[1][1][1], point[2][1][1], point[1][2][1], point[2][2][1] };
    corners[4] = { point[0][0][1], point[1][0][1], point[0][1][1], point[1][1][1],
        point[0][0][2], point[1][0][2], point[0][1][2], point[1][1][2] };
    corners[5] = { point[1][0][1], point[2][0][1], point[1][1][1], point[2][1][1],
        point[1][0][2], point[2][0][2], point[1][1][2], point[2][1][2] };
    corners[6] = { point[0][1][1], point[1][1][1], point[0][2][1], point[1][2][1],
        point[0][1][2], point[1][1][2], point[0][2][2], point[1][2][2] };
    corners[7] = { point[1][1][1], point[2][1][1], point[1][2][1], point[2][2][1],
        point[1][1][2], point[2][1][2], point[1][2][2], point[2][2][2] };
    /*
    std::ofstream cornerCsv("corner.csv");
    for (auto corner : corners) {
        for (int i = 0; i < 8; i++) {
            cornerCsv << corner[i][0] << ",";
            cornerCsv << corner[i][1] << ",";
            cornerCsv << corner[i][2] << "\n";
        }
    }
    cornerCsv.close();
    */

    std::array<int, 9> result = return ninside8Domains(corners, N);
    BOOST_CHECK(result[0] == 0);
    BOOST_CHECK(result[1] == N);
    for (int i = 2; i < 9; i++) {
        BOOST_CHECK(result[i] == 0);
    }
}