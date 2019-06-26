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
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <repa/repa.hpp>

boost::mpi::environment env;

static void test(const boost::mpi::communicator &comm,
                 repa::grids::ParallelLCGrid *grid)
{
    int nlocalcells = grid->n_local_cells();
    BOOST_TEST(nlocalcells >= 0);

    int nglobalcells
        = boost::mpi::all_reduce(comm, nlocalcells, std::plus<int>{});
    BOOST_TEST(nglobalcells > 0);

    auto gs = grid->grid_size();
    BOOST_TEST((nglobalcells == gs[0] * gs[1] * gs[2]));

    BOOST_TEST(grid->n_ghost_cells() <= nglobalcells);
}

BOOST_AUTO_TEST_CASE(test_cell_numbers)
{
    using std::placeholders::_1;
    boost::mpi::communicator comm;
    repa::Vec3d box = {{20., 20., 20.}};
    double mings = 1.0;
    new_test_env(comm, box, mings)
        .with_repart()
        .all_grids()
        .run(std::bind(test, std::cref(comm), _1));
}
