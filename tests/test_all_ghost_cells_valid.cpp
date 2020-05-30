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
 * Tests if all ghost cells have at least one inner neighbor cell,
 * and ensures that each ghost and boundary cells have associated
 * communications.
 */

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ghost_cells_valid
#include <boost/test/unit_test.hpp>

#include "testenv.hpp"
#include <algorithm>
#include <boost/mpi/environment.hpp>
#include <repa/repa.hpp>

using CVBIt = std::vector<bool>::const_iterator;
static bool all_true(const CVBIt &first, const CVBIt &last)
{
    return std::all_of(first, last, [](bool b) { return b; });
}

template <typename Container>
bool contains(Container cont, typename Container::value_type value)
{
    return std::find(std::begin(cont), std::end(cont), value) != std::end(cont);
}

static void test_boundary_has_comm(repa::grids::ParallelLCGrid *grid)
{
    // Checks that boundary cells have an associated sent operation
    auto check_is_send_cell = [](auto c, const auto &gexds) {
        return std::find_if(
                   std::begin(gexds), std::end(gexds),
                   [c](const auto &gexd) { return contains(gexd.send, c); })
               != std::end(gexds);
    };

    auto gexds = grid->get_boundary_info();
    for (const auto c : grid->local_cells()) {
        for (int j = 0; j < 27; ++j) {
            grid->cell_neighbor_index(c, j)
                .visit_if<repa::local_cell_index_type>(
                    [&](repa::local_cell_index_type neighidx) {
                        check_is_send_cell(neighidx, gexds);
                    });
        }
    }
}

static void test_ghost_has_comm(repa::grids::ParallelLCGrid *grid)
{
    // Test that all ghost cells have an associated receive operation
    std::vector<bool> used(grid->ghost_cells().size(), false);

    for (const auto &g : grid->get_boundary_info()) {
        for (const auto &ghost : g.recv) {
            // Ensure valid ghost cell index
            BOOST_TEST(((ghost.value() >= 0)
                        && (ghost.value() < grid->ghost_cells().size())));
            // Each ghost cell can only have one receive operation.
            BOOST_TEST(!used[ghost.value()]);
            used.at(ghost.value()) = true;
        }
    }
    BOOST_TEST(all_true(std::begin(used), std::end(used)));
}

static void test_ghost_has_local(repa::grids::ParallelLCGrid *grid)
{
    // Test that all ghost cells have a neighboring inner cell
    std::vector<bool> used(grid->ghost_cells().size(), false);

    for (const auto c : grid->local_cells()) {
        for (int j = 0; j < 27; ++j) {
            grid->cell_neighbor_index(c, j)
                .visit_if<repa::ghost_cell_index_type>(
                    [&](const auto ghostidx) { used[ghostidx] = true; });
        }
    }
    BOOST_TEST(all_true(std::begin(used), std::end(used)));
}

static void test(const testenv::TEnv &t, repa::grids::ParallelLCGrid *grid)
{
    (void)t;
    test_ghost_has_local(grid);
    test_ghost_has_comm(grid);
    test_boundary_has_comm(grid);
}

BOOST_AUTO_TEST_CASE(test_all_ghost_cells_valid)
{
    testenv::TEnv::default_test_env().with_repart().all_grids().run(test);
}

int main(int argc, char **argv)
{
    boost::mpi::environment mpi_env{argc, argv};
    return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
