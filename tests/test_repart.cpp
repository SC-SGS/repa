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

#define BOOST_TEST_MODULE repart
#include <boost/test/included/unit_test.hpp>

#include "testenv.hpp"
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <repa/repa.hpp>

boost::mpi::environment env;

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

static void test(const TEnv &t, repa::grids::ParallelLCGrid *grid)
{
    (void)t;
    auto nlc = grid->n_local_cells();

    auto all_ones = [nlc]() { return std::vector<double>(nlc, 1.0); };
    auto constant_one = [](int i, int j) { return 1.0; };

    auto cc = CallCounter{};
    bool grid_changed;
    BOOST_CHECK_NO_THROW(
        grid_changed = grid->repartition(all_ones, constant_one, std::ref(cc)));

    // Callback needs to be evaluated exactly once and only if grid was changed.
    BOOST_TEST(if_then(grid_changed, cc.get_count() == 1));
    BOOST_TEST(if_then(!grid_changed, cc.get_count() == 0));
}

BOOST_AUTO_TEST_CASE(test_repart)
{
    boost::mpi::environment env;
    default_test_env().without_repart().all_grids().run(test);
}
