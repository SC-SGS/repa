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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE repart
#include <boost/test/unit_test.hpp>

#include "testenv.hpp"
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <repa/repa.hpp>

struct CallCounter {
    CallCounter() : cnt(0)
    {
    }
    void operator()()
    {
        cnt++;
    }
    int get_count()
    {
        return cnt;
    }

private:
    int cnt;
};

bool if_then(bool a, bool b)
{
    return !a || b;
}

static void test(const testenv::TEnv &t, repa::grids::ParallelLCGrid *grid)
{
    (void)t;
    auto f_all_ones = [](size_t size) {
        return [size]() { return std::vector<double>(size, 1.0); };
    };
    auto constant_one
        = [](repa::local_cell_index_type i,
             repa::local_or_ghost_cell_index_type j) { return 1.0; };

    auto cc = CallCounter{};
    bool grid_changed = false;
    BOOST_CHECK_NO_THROW(
        grid_changed = grid->repartition(f_all_ones(grid->local_cells().size()),
                                         constant_one, std::ref(cc)));

    // Callback needs to be evaluated exactly once and only if grid was changed.
    BOOST_TEST(if_then(grid_changed, cc.get_count() == 1));
    BOOST_TEST(if_then(!grid_changed, cc.get_count() == 0));
}

BOOST_AUTO_TEST_CASE(test_repart)
{
    testenv::TEnv::default_test_env().without_repart().all_grids().run(test);
}

int main(int argc, char **argv)
{
    boost::mpi::environment mpi_env{argc, argv};
    return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
