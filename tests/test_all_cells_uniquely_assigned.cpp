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
#define BOOST_TEST_MODULE all_cells_assigned
#include <boost/test/unit_test.hpp>

#include "testenv.hpp"
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <random>
#include <repa/repa.hpp>

static bool is_assigned(const boost::mpi::communicator &comm,
                        repa::grids::ParallelLCGrid *grid,
                        const repa::Vec3d &pos)
{
    repa::local_cell_index_type cellidx;
    bool has_cell = false;
    try {
        cellidx = grid->position_to_cell_index(pos);
        has_cell = true;
    }
    catch (...) {
    }

    // Ensure a process claiming to be responsibe is indeed the owner.
    if (has_cell) {
        BOOST_CHECK(grid->position_to_rank(pos) == comm.rank());
        BOOST_CHECK(
            (cellidx >= 0
             && static_cast<size_t>(cellidx) < grid->local_cells().size()));
    }
    else {
        try {
            // Might throw which is okay. But if it returns, it must return
            // the rank of the calling process.
            auto owner = grid->position_to_rank(pos);
            BOOST_CHECK(owner != comm.rank());
        }
        catch (...) {
        }
    }

    return has_cell;
}

static void test(const testenv::TEnv &t, repa::grids::ParallelLCGrid *grid)
{
    // Test reguar grid
    auto cell_size = grid->cell_size();
    repa::Vec3d pos;
    std::vector<int> nowners;
    nowners.reserve(grid->grid_size()[0] * grid->grid_size()[1]
                    * grid->grid_size()[2]);

    for (pos[0] = .5 * cell_size[0]; pos[0] < t.box()[0];
         pos[0] += cell_size[0]) {
        for (pos[1] = .5 * cell_size[1]; pos[1] < t.box()[1];
             pos[1] += cell_size[1]) {
            for (pos[2] = .5 * cell_size[2]; pos[2] < t.box()[2];
                 pos[2] += cell_size[2]) {
                nowners.push_back(is_assigned(t.comm(), grid, pos) ? 1 : 0);
            }
        }
    }

    std::vector<int> nowners_global;
    nowners_global.reserve(nowners.size());
    boost::mpi::reduce(t.comm(), nowners, nowners_global, std::plus{}, 0);
    if (t.comm().rank() == 0) {
        BOOST_CHECK(std::all_of(nowners_global.begin(), nowners_global.end(),
                                [](int i) { return i == 1; }));
    }
}

// Gridbased and Diffusion do not allow for position_to_rank after
// repartitioning, so test statically.
BOOST_AUTO_TEST_CASE(test_all_cells_uniquely_assigned)
{
    testenv::TEnv::default_test_env().with_repart().all_grids().run(test);
}

int main(int argc, char **argv)
{
    boost::mpi::environment mpi_env{argc, argv};
    return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
