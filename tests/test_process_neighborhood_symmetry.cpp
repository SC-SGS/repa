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
 * Checks the symmetry of the process neighborhood relation.
 */

#define BOOST_TEST_MODULE process_neighborhood_symmetry
#include <boost/test/included/unit_test.hpp>

#include "testenv.hpp"
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/vector.hpp>
#include <repa/repa.hpp>

static void test(const testenv::TEnv &t, repa::grids::ParallelLCGrid *grid)
{
    const auto &comm = t.comm();
    std::vector<int> neighranks;
    neighranks.reserve(grid->n_neighbors());

    // Test validity and uniqueness (correct size).
    for (repa::rank_index_type i = 0; i < grid->n_neighbors(); ++i) {
        auto n = grid->neighbor_rank(i);

        BOOST_TEST(((n >= 0) && (n < comm.size())));

        for (const auto n0 : neighranks) {
            BOOST_TEST(n != n0);
        }
        neighranks.push_back(n);
    }

    std::vector<decltype(neighranks)> nrankss;
    boost::mpi::all_gather(comm, neighranks, nrankss);

    // Check symmetry of process neighborhoods
    for (int rank1 = 0; rank1 < comm.size(); ++rank1) {
        for (auto rank2 : nrankss[rank1]) {
            auto it = std::find(std::begin(nrankss[rank2]),
                                std::end(nrankss[rank2]), rank1);
            BOOST_TEST((it != std::end(nrankss[rank2])));
        }
    }
}

BOOST_AUTO_TEST_CASE(test_process_neighborhood_symmetry)
{
    boost::mpi::environment env;
    testenv::TEnv::default_test_env().with_repart().all_grids().run(test);
}
