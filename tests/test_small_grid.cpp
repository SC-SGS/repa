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
 * Tests if the grids can handle a small grid with dimensions <= 1.0.
 * Mainly for gridbased to test tetra::precision.
 */

#define BOOST_TEST_MODULE small_grid
#include <boost/test/included/unit_test.hpp>

#include "testenv.hpp"
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <random>
#include <repa/repa.hpp>

static void test(const testenv::TEnv &t, repa::grids::ParallelLCGrid *grid)
{
    BOOST_TEST(grid->n_local_cells() > 0);
}

// Gridbased and Diffusion do not allow for position_to_rank after
// repartitioning, so test statically.
BOOST_AUTO_TEST_CASE(test_small_grid)
{
    boost::mpi::environment env;
    testenv::TEnv::custom_test_env({1.0, 1.0, 1.0}, 0.162)
        .without_repart()
        .all_grids()
        .run(test);
}
