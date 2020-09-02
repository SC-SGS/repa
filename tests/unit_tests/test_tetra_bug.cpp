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

/** This test case tests for several scenarios that revealed past tetra errors.
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tetra_bug
#include <boost/test/unit_test.hpp>

#include <array>

#include <repa/grids/util/tetra.hpp>
#include <repa/grids/util/vec_arith.hpp>

struct TestCase {
    int nproc;
    repa::Vec3i dims;
    repa::Vec3d box_size;
    double min_gs;
    repa::Vec3i grid_size;
    repa::Vec3d cell_size;
    std::vector<repa::Vec3d> gridpoints;
};

/** Gridpoints are saved in OpenMPI linearization
 * (MPI_Cart_rank and MPI_Cart_coord).
 * The following two functions resemble them.
 */
static repa::Vec3i unlinearize(int r, const repa::Vec3i &dims)
{
    repa::Vec3i res;
    res[2] = r % dims[2];
    r /= dims[2];
    res[1] = r % dims[1];
    r /= dims[1];
    assert(r < dims[0]);
    res[0] = r;
    return res;
}

static int linearize(repa::Vec3i r, const repa::Vec3i &dims)
{
    return r[0] * dims[1] * dims[2] + r[1] * dims[2] + r[2];
}

/** Copy of gridbased.cpp: GridBasedGrid::unshifted_bounding_box(rank_type).
 */
repa::util::tetra::BoundingBox
bounding_box(int which_subdomain,
             const std::vector<repa::Vec3d> &gridpoints,
             const repa::Vec3i &dims)
{
    const auto coord = unlinearize(which_subdomain, dims);

    std::array<repa::Vec3d, 8> ps;
    std::array<repa::Vec3i, 8> ms;
    size_t i = 0;
    // Ranks holding the bounding box grid points of "r" = (c0, c1, c2) are:
    // (c0,     c1,     c2) upper right back corner,
    // (c0 - 1, c1,     c2) upper left back corner,
    // (c0,     c1 - 1, c2) lower right back corner,
    // (c0 - 1, c1 - 1, c2) lower left back corner
    // (c0,     c1,     c2 - 1) upper right front corner,
    // ... 2 more ...
    // (c0 - 1, c1 - 1, c2 - 1) lower left front corner
    // In total the set: {c0, c0 - 1} x {c1, c1 - 1} x {c2, c2 - 1}
    repa::Vec3i off;
    for (off[2] = 0; off[2] <= 1; ++off[2]) {
        for (off[1] = 0; off[1] <= 1; ++off[1]) {
            for (off[0] = 0; off[0] <= 1; ++off[0]) {
                using namespace repa::util::vector_arithmetic;
                repa::Vec3i nc = (coord - off) % dims;
                int proc = linearize(nc, dims);

                // Mirror the gridpoint back to where this subdomain is
                // expecting it.
                const repa::Vec3i mirror
                    = -static_cast_vec<repa::Vec3i>((coord == 0) && (off == 1));
                ps[i] = gridpoints[proc];
                ms[i] = mirror;
                i++;
            }
        }
    }
    return repa::util::tetra::BoundingBox{std::move(ps), std::move(ms)};
}

void test_impl(const TestCase &tc)
{
    repa::util::tetra::init_tetra(tc.min_gs, tc.box_size);
    using namespace repa::util::tetra;

    std::vector<Octagon> octas;

    for (int i = 0; i < tc.nproc; ++i) {
        octas.emplace_back(bounding_box(i, tc.gridpoints, tc.dims), tc.min_gs);
        BOOST_TEST(octas[i].is_valid());
    }

    // Iterate over all centers of cells of width = min_gs x min_gs x min_gs.
    bool all_accepted_exactly_once = true;
    repa::Vec3d p;
    for (p[0] = tc.min_gs / 2; p[0] < tc.box_size[0]; p[0] += tc.min_gs) {
        for (p[1] = tc.min_gs / 2; p[1] < tc.box_size[1]; p[1] += tc.min_gs) {
            for (p[2] = tc.min_gs / 2; p[2] < tc.box_size[2];
                 p[2] += tc.min_gs) {

                const auto naccept = std::count_if(
                    octas.begin(), octas.end(),
                    [p](const Octagon &o) { return o.contains(p); });
                if (naccept != 1) {
                    std::cerr << p[0] << ", " << p[1] << ", " << p[2] << ": "
                              << naccept << std::endl;
                    all_accepted_exactly_once = false;
                }
            }
        }
    }
    BOOST_TEST(all_accepted_exactly_once);
}

#define REGISTER_TETRA_TEST(testcase)                                          \
    BOOST_AUTO_TEST_CASE(test_tetra_bug_##testcase)                            \
    {                                                                          \
        test_impl(testcase);                                                   \
    }

// Again, C++ prior to C++20 does not support designated initializers, which
// would be a massive increase in readability...
// We add them anyway but define them away for now.
// NOTE that this REQUIRES the initializers to be in the exact order of the
// declaration in struct TestCase.
// And yes, __ is a reserved identifier...
#define __(x)

/** Checks if tetra can handle a certain subdomain layout.
 * This subdomain layout has gridpoints shifted over the periodic boundary.
 * It is a case which is known to fail before commit ce1d346.
 */
static const TestCase test_case_pre_ce1d346 = {
    __(.nproc =) 8,
    __(.dims =) repa::Vec3i{2, 2, 2},
    __(.box_size =) repa::Vec3d{80., 80., 80.},
    __(.min_gs =) 1.,
    __(.grid_size =) repa::Vec3i{80, 80, 80},
    __(.cell_size =) repa::Vec3d{1., 1., 1.},
    __(.gridpoints =){{{35.567, 38.4418, 45.0975}},
                      {{35.7365, 38.3409, 75.8135}},
                      {{35.9562, 82.1775, 44.5612}},
                      {{35.9879, 82.2762, 76.3616}},
                      {{83.8035, 38.3232, 45.2142}},
                      {{83.641, 38.3249, 75.6815}},
                      {{83.4199, 82.2701, 44.6563}},
                      {{83.3893, 82.2873, 76.2623}}},
};

REGISTER_TETRA_TEST(test_case_pre_ce1d346)

/** This subdomain layout has gridpoints shifted over the periodic boundary so
 * slightly that in conjunction with the specific box_size causes rounding
 * errors if vertices, etc. are shifted as doubles.
 * It is a case which is known to fail before commit ce1d346.
 */
static const TestCase test_case_pre_183a75b = {
    __(.nproc =) 24,
    __(.dims =) repa::Vec3i{4, 3, 2},
    __(.box_size =)
        repa::Vec3d{133.88659001643387, 133.88659001643387, 133.88659001643387},
    __(.min_gs =) 2.5,
    __(.grid_size =) repa::Vec3i{53, 53, 53},
    __(.cell_size =) repa::Vec3d{2.52616207578177113207, 2.52616207578177113207,
                                 2.52616207578177113207},
    __(.gridpoints
       =){{{43.494856734183919, 53.899425582781973, 67.98611172561732}},
          {{38.971624032523593, 47.820236497391079, 133.64229606120034}},
          {{48.242421542159377, 77.768214283173549, 68.473723891871259}},
          {{38.786088957104226, 85.201562482732029, 133.45066967165698}},
          {{33.471647504108468, 133.88659001643387, 66.943295008216936}},
          {{33.471647504108468, 133.88659001643387, 133.88659001643387}},
          {{67.250249064965928, 59.973796588371187, 67.172297899238657}},
          {{66.898486932319869, 47.721358323519432, 134.31461522354823}},
          {{66.405435009323256, 77.51303233813988, 65.731781630025282}},
          {{66.875322544805471, 86.07079483788695, 134.99311345583953}},
          {{66.943295008216936, 133.88659001643387, 66.943295008216936}},
          {{66.943295008216936, 133.88659001643387, 133.88659001643387}},
          {{88.646807832764424, 59.716307749073508, 66.626508470926765}},
          {{95.318981994171224, 48.191191477600498, 134.61051882204632}},
          {{86.238282618015134, 75.421278692274157, 67.74850013025619}},
          {{94.887030555244323, 85.372682118099746, 135.25482122377537}},
          {{100.41494251232541, 133.88659001643387, 66.943295008216936}},
          {{100.41494251232541, 133.88659001643387, 133.88659001643387}},
          {{133.72006927060264, 44.708757714230536, 66.337285296722158}},
          {{133.69913693065274, 44.704261775564298, 134.48850223536112}},
          {{133.65458952221903, 89.346596944554861, 66.85129813517382}},
          {{133.66110018751812, 89.331147167707869, 133.9729645568525}},
          {{133.88659001643387, 133.88659001643387, 66.943295008216936}},
          {{133.88659001643387, 133.88659001643387, 133.88659001643387}}},
};

REGISTER_TETRA_TEST(test_case_pre_183a75b)
