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
 * Checks the full shell neighborhood symmetry for inner cells.
 */

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE neighborhood_symmetry
#include <boost/test/unit_test.hpp>

#include "testenv.hpp"
#include <boost/mpi/environment.hpp>
#include <repa/repa.hpp>

static bool has_neighbor(repa::grids::ParallelLCGrid *grid,
                         repa::local_cell_index_type d,
                         repa::local_cell_index_type c)
{
    for (int j = 0; j < 27; ++j) {
        if (grid->cell_neighbor_index(d, j) == c)
            return true;
    }
    return false;
}

static void test(const testenv::TEnv &t, repa::grids::ParallelLCGrid *grid)
{
    (void)t;

    const auto nlocalcells = grid->local_cells().size();
    const auto nghostcells = grid->ghost_cells().size();

    for (const auto c : grid->local_cells()) {
        for (int j = 0; j < 27; ++j) {
            const auto d = grid->cell_neighbor_index(c, j);

            // If "d" is a inner cell and neighbors "c", then "c" must also
            // neighbor "d".
            if (d.is<repa::local_cell_index_type>()) {
                const auto li = d.as<repa::local_cell_index_type>();
                BOOST_CHECK(li >= 0);
                BOOST_CHECK(static_cast<size_t>(li) < nlocalcells);
                BOOST_CHECK(
                    li < static_cast<repa::local_cell_index_type::value_type>(
                        nlocalcells));
                BOOST_CHECK(
                    has_neighbor(grid, d.as<repa::local_cell_index_type>(), c));
            }
            else {
                const auto gi = d.as<repa::ghost_cell_index_type>();
                BOOST_CHECK(gi >= 0);
                BOOST_CHECK(static_cast<size_t>(gi) < nghostcells);
                BOOST_CHECK(
                    gi < static_cast<repa::ghost_cell_index_type::value_type>(
                        nghostcells));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_neighborhood_symmetry)
{
    testenv::TEnv::default_test_env().with_repart().all_grids().run(test);
}

int main(int argc, char **argv)
{
    boost::mpi::environment mpi_env{argc, argv};
    return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
