/**
 * Copyright 2020 The repa authors
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
 * Checks if tetra can handle a certain subdomain layout.
 * This subdomain layout has gridpoints shifted over the periodic boundary.
 * It is a case which is known to fail before commit ce1d346.
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE tetra
#include <boost/test/unit_test.hpp>

#include <array>

#include <repa/grids/util/tetra.hpp>
#include <repa/grids/util/vec_arith.hpp>

/** Data of erroneous run.
 */
static const int nproc = 8;
static const auto dims = repa::Vec3i{2, 2, 2};
static const auto box_size = repa::Vec3d{80., 80., 80.};
static const double min_gs = 1.;
static const auto grid_size = repa::Vec3i{80, 80, 80};
static const auto cell_size = repa::Vec3d{1., 1., 1.};
static const std::array<repa::Vec3d, nproc> gridpoints
    = {{{{35.567, 38.4418, 45.0975}},
        {{35.7365, 38.3409, 75.8135}},
        {{35.9562, 82.1775, 44.5612}},
        {{35.9879, 82.2762, 76.3616}},
        {{83.8035, 38.3232, 45.2142}},
        {{83.641, 38.3249, 75.6815}},
        {{83.4199, 82.2701, 44.6563}},
        {{83.3893, 82.2873, 76.2623}}}};

/** Gridpoints are saved in OpenMPI linearization
 * (MPI_Cart_rank and MPI_Cart_coord).
 * The following two functions resemble them.
 */
static repa::Vec3i unlinearize(int r)
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

static int linearize(repa::Vec3i r)
{
    return r[0] * dims[1] * dims[2] + r[1] * dims[2] + r[2];
}

/** Copy of gridbased.cpp: GridBasedGrid::unshifted_bounding_box(rank_type).
 */
repa::util::tetra::BoundingBox bounding_box(int which_subdomain)
{
    const auto coord = unlinearize(which_subdomain);

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
                int proc = linearize(nc);

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

BOOST_AUTO_TEST_CASE(test_tetra_bug)
{
    repa::util::tetra::init_tetra(min_gs, box_size);
    using namespace repa::util::tetra;

    std::array<Octagon, nproc> octas;

    for (int i = 0; i < nproc; ++i) {
        octas[i] = Octagon(bounding_box(i), min_gs);
        BOOST_TEST(octas[i].is_valid());
    }

    // Iterate over all centers of cells of width = min_gs x min_gs x min_gs.
    bool all_accepted_exactly_once = true;
    repa::Vec3d p;
    for (p[0] = min_gs / 2; p[0] < box_size[0]; p[0] += min_gs) {
        for (p[1] = min_gs / 2; p[1] < box_size[1]; p[1] += min_gs) {
            for (p[2] = min_gs / 2; p[2] < box_size[2]; p[2] += min_gs) {

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
