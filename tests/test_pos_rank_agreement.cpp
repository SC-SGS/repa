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

#define BOOST_TEST_MODULE pos_rank_agreement
#include <boost/test/included/unit_test.hpp>

#include "testenv.hpp"
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/vector.hpp>
#include <random>
#include <repa/repa.hpp>

template <typename T>
static void test_agreement(const boost::mpi::communicator &comm, T value)
{
    std::vector<T> values;
    boost::mpi::all_gather(comm, value, values);

    for (const auto v : values) {
        BOOST_TEST(v == value);
    }
}

template <typename E>
bool ignore_message(const E &e)
{
    (void)e;
    return true;
}

static std::vector<int> neighranks(repa::grids::ParallelLCGrid *grid)
{
    std::vector<int> res;
    res.reserve(grid->n_neighbors());
    for (int i = 0; i < grid->n_neighbors(); ++i)
        res.push_back(grid->neighbor_rank(i));
    return res;
}

static void test_position(const boost::mpi::communicator &comm,
                          repa::grids::ParallelLCGrid *grid,
                          const repa::Vec3d &pos)
{
    int rank = grid->position_to_rank(pos.data());

    // Check correctness of "rank"
    BOOST_TEST(((0 <= rank) && (rank < comm.size())));

    // Additional checks for owner of "pos"
    if (rank == comm.rank()) {
        // We cannot check that grid->position_to_neighidx(pos.data())
        // does not throw, here. A grid could implement a "full halo" as ghost
        // layer and need self-communication (e.g. cart grid with any direction
        // only having a node_grid of 1).

        // Check validity of resolved cell
        int cidx = -1;
        BOOST_CHECK_NO_THROW(cidx = grid->position_to_cell_index(pos.data()));

        BOOST_TEST(((0 <= cidx) && (cidx < grid->n_local_cells())));
    }
    else {
        // Check that all other processes refuse to resolve "pos"
        // to a cell index.
        BOOST_CHECK_EXCEPTION(grid->position_to_cell_index(pos.data()),
                              std::domain_error, ignore_message);
    }

    // Check that position_to_neighidx can resolve "pos" if it actually is on a
    // neighboring process
    auto neighborranks
        = neighranks(grid); // Okay, recreating this vecor for every position is
                            // stupid, but I currently don't care.
    if (std::find(std::begin(neighborranks), std::end(neighborranks), rank)
        != std::end(neighborranks)) {
        int nidx = -1;
        BOOST_CHECK_NO_THROW(nidx = grid->position_to_neighidx(pos.data()));
        BOOST_TEST(((nidx >= 0) && (nidx <= grid->n_neighbors())));
    }

    test_agreement(comm, rank);
}

static void test(const TEnv &t, repa::grids::ParallelLCGrid *grid)
{
    const auto &comm = t.comm;
    const auto &box = t.box;
    // Test reguar grid
    auto cell_size = grid->cell_size();
    repa::Vec3d pos;

    for (pos[0] = 0.1; pos[0] < box[0]; pos[0] += .9 * cell_size[0]) {
        for (pos[1] = 0.2; pos[1] < box[1]; pos[1] += .7 * cell_size[1]) {
            for (pos[2] = 0.3; pos[2] < box[2]; pos[2] += 1.3 * cell_size[2]) {
                test_position(comm, grid, pos);
            }
        }
    }

    // Test some random positions
    std::random_device rd;
    std::mt19937 gen(rd());
    typedef std::uniform_real_distribution<> Dist;
    std::array<Dist, 3> diss{Dist{0.0, box[0]}, Dist{0.0, box[1]},
                             Dist{0.0, box[2]}};

    for (size_t i = 0; i < 10'000; ++i) {
        if (comm.rank() == 0) {
            for (size_t d = 0; d < pos.size(); ++d)
                pos[d] = diss[d](gen);
        }
        boost::mpi::broadcast(comm, pos, 0);

        test_position(comm, grid, pos);
    }
}

BOOST_AUTO_TEST_CASE(test_pos_rank_agreement)
{
    boost::mpi::environment env;
    // Gridbased and Diffusion do not allow for position_to_rank after
    // repartitioning, so test statically.
    default_test_env()
        .without_repart()
        .only({repa::GridType::DIFF, repa::GridType::GRIDBASED,
               repa::GridType::HYB_GP_DIFF})
        .run(test);
    default_test_env()
        .with_repart()
        .all_grids()
        .exclude({repa::GridType::DIFF, repa::GridType::GRIDBASED,
                  repa::GridType::HYB_GP_DIFF})
        .run(test);
}