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
 * Checks for the correct creation of supported pargrids.
 */

#define BOOST_TEST_MODULE hybrid_method_switch
#include <boost/test/included/unit_test.hpp>

#include "testenv.hpp"
#include <boost/mpi/environment.hpp>
#include <repa/repa.hpp>

static void set_toggle(const testenv::TEnv &t,
                       repa::grids::ParallelLCGrid *grid)
{
    (void)t;
    grid->command("toggle");
}

static void set_set(const testenv::TEnv &t, repa::grids::ParallelLCGrid *grid)
{
    static int ncalls = 0;
    (void)t;

    if (ncalls % 2 == 0)
        grid->command("set diffusion");
    else
        grid->command("set graph");
    ncalls++;
}

BOOST_AUTO_TEST_CASE(test_hybrid_method_switch)
{
    boost::mpi::environment env;
    testenv::TEnv::default_test_env()
        .with_repart_twice()
        .only({repa::GridType::HYB_GP_DIFF})
        .run(set_toggle);

    testenv::TEnv::default_test_env()
        .with_repart_twice()
        .only({repa::GridType::HYB_GP_DIFF})
        .run(set_set);
}
