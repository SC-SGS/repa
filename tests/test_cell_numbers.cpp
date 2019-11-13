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
 * Checks if the number of local cells on each process is meaningful.
 */

#define BOOST_TEST_MODULE cell_numbers
#include <boost/test/included/unit_test.hpp>

#include "testenv.hpp"
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/environment.hpp>
#include <repa/repa.hpp>


static bool if_then(bool a, bool b)
{
    return !a || b;
}

static void test(const testenv::TEnv &t, repa::grids::ParallelLCGrid *grid)
{
    int nlocalcells = grid->n_local_cells();
    BOOST_TEST(nlocalcells >= 0);

    int nglobalcells
        = boost::mpi::all_reduce(t.comm, nlocalcells, std::plus<int>{});
    BOOST_TEST(nglobalcells > 0);

    auto gs = grid->grid_size();
    BOOST_TEST((nglobalcells == gs[0] * gs[1] * gs[2]));

    // Full-halo grids might have more ghost than local cells on 1 process only
    // and very small grids.
    BOOST_TEST(
        if_then(t.comm.size() >= 2, grid->n_ghost_cells() <= nglobalcells));
}

BOOST_AUTO_TEST_CASE(test_cell_numbers)
{
    boost::mpi::environment env;
    default_test_env().with_repart().all_grids().run(test);
}
