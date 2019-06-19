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
 * Tests if all ghost cells have at least one inner neighbor cell.
 */
#define BOOST_TEST_MODULE ghost_cells_valid
#include <boost/test/included/unit_test.hpp>

#include <algorithm>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <repa/repa.hpp>

static void test(repa::grids::ParallelLCGrid *grid)
{
    std::vector<bool> used(grid->n_ghost_cells(), false);

    for (int c = 0; c < grid->n_local_cells(); ++c) {
        for (int j = 0; j < 27; ++j) {
            int d = grid->cell_neighbor_index(c, j);

            if (d >= grid->n_local_cells()) {
                used[d - grid->n_local_cells()] = true;
            }
        }
    }

    auto all_true = [](auto first, auto last) {
        return std::all_of(first, last, [](bool b) { return b; });
    };
    BOOST_TEST(all_true(std::begin(used), std::end(used)));
}

BOOST_AUTO_TEST_CASE(test_all_ghost_cells_valid)
{
    boost::mpi::environment env;
    boost::mpi::communicator comm;

    for (const auto gt : repa::supported_grid_types()) {
        if (comm.rank() == 0) {
            std::cout << "Checking grid '" << repa::grid_type_to_string(gt)
                      << "'" << std::endl;
        }
        auto up = repa::make_pargrid(gt, comm, {{20., 20., 20.}}, 1.0);
        test(up.get());
    }
}
