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
 * Tests if all processes agree on pos_to_rank
 */

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE extra_params
#include <boost/test/unit_test.hpp>

#include "testenv.hpp"
#include <boost/mpi/environment.hpp>
#include <repa/repa.hpp>

/** Helper struct that counts the number of times ep.subdomain_center
 * has been called.
 */
static struct EPCallCount {
    int ncalls = 0;
    repa::ExtraParams ep;

    std::pair<int, repa::Vec3d> operator()(repa::local_cell_index_type)
    {
        ncalls++;
        return std::make_pair(0, repa::Vec3d{0., 0., 0.});
    }

    EPCallCount()
        : ep(repa::ExtraParams{
            std::bind(&EPCallCount::operator(), this, std::placeholders::_1),
            {}})
    {
    }
} epcallcount;

static void test(const testenv::TEnv &t, repa::grids::ParallelLCGrid *grid)
{
    static int ncalls = 0;
    // On second call (the call to 'test' after repartition),
    // ExtraParams::subdomain_midpoint must have been called for each cell
    // exactly once. Before that it must not have been called.
    BOOST_CHECK(epcallcount.ncalls == ncalls * grid->n_local_cells());
    ncalls++;
}

// Only check for gridbased since it is currently the only grid using
// ExtraParams
BOOST_AUTO_TEST_CASE(test_extra_params)
{
    testenv::TEnv::default_test_env(epcallcount.ep)
        .with_repart()
        .only({repa::GridType::GRIDBASED})
        .run(test);
}

int main(int argc, char **argv)
{
    boost::mpi::environment mpi_env{argc, argv};
    return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
