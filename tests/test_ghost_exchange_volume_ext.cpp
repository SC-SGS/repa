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
 * Tests the symmetry of ghost exchange across processes.
 */
#define BOOST_TEST_MODULE ghost_exchange_volume
#include <boost/test/included/unit_test.hpp>

#include "testenv.hpp"
#include <algorithm>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/vector.hpp>
#include <repa/repa.hpp>


// Serialization for GhostExchangeDesc in order to gather and check them.
namespace boost {
namespace serialization {
template <typename Archive>
void load(Archive &ar,
          repa::grids::GhostExchangeDesc &g,
          const unsigned int /* file_version */)
{
    ar >> g.dest;
    ar >> g.recv;
    ar >> g.send;
}

template <typename Archive>
void save(Archive &ar,
          const repa::grids::GhostExchangeDesc &g,
          const unsigned int /* file_version */)
{
    ar << g.dest;
    ar << g.recv;
    ar << g.send;
}

template <class Archive>
void serialize(Archive &ar,
               repa::grids::GhostExchangeDesc &g,
               const unsigned int file_version)
{
    split_free(ar, g, file_version);
}
} // namespace serialization
} // namespace boost

static void test(const TEnv &t, repa::grids::ParallelLCGrid *grid)
{
    const auto &comm = t.comm;
    auto gexds = grid->get_boundary_info();
    auto idx_to_glo = [grid](int idx) { return grid->global_hash(idx); };

    // Transform to global hashes
    for (auto &g : gexds) {
        std::transform(std::begin(g.recv), std::end(g.recv), std::begin(g.recv),
                       idx_to_glo);
        std::transform(std::begin(g.send), std::end(g.send), std::begin(g.send),
                       idx_to_glo);
    }

    std::vector<decltype(gexds)> gexdss;
    boost::mpi::all_gather(comm, gexds, gexdss);

    auto find_comm = [](const std::vector<repa::grids::GhostExchangeDesc> &gs,
                        int rank) -> const repa::grids::GhostExchangeDesc& {
        auto it
            = std::find_if(std::begin(gs), std::end(gs),
                           [rank](const auto &g) { return g.dest == rank; });
        BOOST_TEST((it != std::end(gs)));
        return *it;
    };

    // Check if send indices fit receive indices on the other side.
    for (int r = 0; r < comm.size(); ++r) {
        for (const auto &rg : gexdss[r]) {
            const auto &counterpart = find_comm(gexdss[rg.dest], r);

            for (size_t i = 0; i < rg.send.size(); ++i) {
                auto sglo = rg.send[i];
                auto rglo = counterpart.recv[i];

                BOOST_TEST(sglo == rglo);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_ghost_exchange_volume)
{
    boost::mpi::environment env;
    default_test_env().without_repart().all_grids().run(test);
}
