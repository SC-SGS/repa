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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE ghost_exchange_volume
#include <boost/test/unit_test.hpp>

#include "testenv.hpp"
#include <algorithm>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/serialization/vector.hpp>
#include <repa/repa.hpp>

/** Analogously to repa::GhostExchangeDecv, however, stores all indices as
 * global ones in order to compare them across processes.
 */
struct GlobalizedGhostExchangeDesc {
    repa::rank_type dest;
    std::vector<repa::global_cell_index_type> recv;
    std::vector<repa::global_cell_index_type> send;

    GlobalizedGhostExchangeDesc() : dest(-1)
    {
    }
};

// Serialization for GlobalizedGhostExchangeDesc in order to gather and check
// them.
namespace boost {
namespace serialization {
template <typename Archive>
void load(Archive &ar,
          GlobalizedGhostExchangeDesc &g,
          const unsigned int /* file_version */)
{
    ar >> g.dest;
    ar >> g.recv;
    ar >> g.send;
}

template <typename Archive>
void save(Archive &ar,
          const GlobalizedGhostExchangeDesc &g,
          const unsigned int /* file_version */)
{
    ar << g.dest;
    ar << g.recv;
    ar << g.send;
}

template <class Archive>
void serialize(Archive &ar,
               GlobalizedGhostExchangeDesc &g,
               const unsigned int file_version)
{
    split_free(ar, g, file_version);
}
} // namespace serialization
} // namespace boost

static void test(const testenv::TEnv &t, repa::grids::ParallelLCGrid *grid)
{
    const auto &comm = t.comm();
    auto idx_to_glo = [&grid](const auto idx) {
        return grid->global_hash(repa::local_or_ghost_cell_index_type{idx});
    };
    std::vector<GlobalizedGhostExchangeDesc> ggexds;

    {
        const auto gexds = grid->get_boundary_info();
        ggexds.resize(gexds.size());

        // Transform to global hashes
        for (size_t i = 0; i < gexds.size(); ++i) {
            ggexds[i].dest = gexds[i].dest;
            ggexds[i].recv.reserve(gexds[i].recv.size());
            ggexds[i].send.reserve(gexds[i].send.size());
            std::transform(std::begin(gexds[i].recv), std::end(gexds[i].recv),
                           std::back_inserter(ggexds[i].recv), idx_to_glo);
            std::transform(std::begin(gexds[i].send), std::end(gexds[i].send),
                           std::begin(ggexds[i].send), idx_to_glo);
        }
    }

    std::vector<decltype(ggexds)> ggexdss;
    boost::mpi::all_gather(comm, ggexds, ggexdss);

    auto find_comm = [](const std::vector<GlobalizedGhostExchangeDesc> &gs,
                        int rank) -> const GlobalizedGhostExchangeDesc & {
        auto it
            = std::find_if(std::begin(gs), std::end(gs),
                           [rank](const auto &g) { return g.dest == rank; });
        BOOST_TEST((it != std::end(gs)));
        return *it;
    };

    // Check if send indices fit receive indices on the other side.
    for (int r = 0; r < comm.size(); ++r) {
        for (const auto &rg : ggexdss[r]) {
            const auto &counterpart = find_comm(ggexdss[rg.dest], r);

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
    testenv::TEnv::default_test_env().without_repart().all_grids().run(test);
}

int main(int argc, char **argv)
{
    boost::mpi::environment mpi_env{argc, argv};
    return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
