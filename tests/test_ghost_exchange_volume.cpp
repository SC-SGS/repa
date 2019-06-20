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

static bool if_then(bool b1, bool b2)
{
    return b2 || !b1;
}

static std::vector<int> neighranks(repa::grids::ParallelLCGrid *grid)
{
    std::vector<int> res;
    res.reserve(grid->n_neighbors());
    for (int i = 0; i < grid->n_neighbors(); ++i)
        res.push_back(grid->neighbor_rank(i));
    return res;
}

static void test(const boost::mpi::communicator &comm,
                 repa::grids::ParallelLCGrid *grid)
{
    auto gexds = grid->get_boundary_info();
    auto neighborranks = neighranks(grid);

    // Verify consistency of neighbor information with ghost communications
    for (auto rank : neighborranks) {
        BOOST_TEST(
            (std::find_if(std::begin(gexds), std::end(gexds),
                          [rank](auto gexd) { return gexd.dest == rank; })
             != std::end(gexds)));
    }

    // Validity of exchange descriptors
    for (const auto &g : gexds) {
        BOOST_TEST((g.dest >= 0 && g.dest < comm.size()));
        BOOST_TEST(g.recv.size() > 0);
        BOOST_TEST(g.send.size() > 0);
    }

    // Validity of cell indices
    for (const auto &g : gexds) {
        // Verify consistency of ghost communications with neighbors
        BOOST_TEST((std::find(std::begin(neighborranks),
                              std::end(neighborranks), g.dest)
                    != std::end(neighborranks)));
        for (auto sendc : g.send) {
            BOOST_TEST(((0 <= sendc) && (sendc < grid->n_local_cells())));
        }
        for (auto recvc : g.recv) {
             BOOST_TEST(((grid->n_local_cells() <= recvc)
                        && (recvc < grid->n_local_cells()
                                        + grid->n_ghost_cells())));
        }
    }

    // Note, although gathered on all processes, all indices will still be in
    // local to the respective process they came from.
    std::vector<decltype(gexds)> gexdss;
    boost::mpi::all_gather(comm, gexds, gexdss);

    auto find_comm = [](std::vector<repa::grids::GhostExchangeDesc> &gs,
                        int rank) {
        auto it
            = std::find_if(std::begin(gs), std::end(gs),
                           [rank](const auto &g) { return g.dest == rank; });
        BOOST_TEST((it != std::end(gs)));
        return *it;
    };

    // Check if send indices fit receive indices on the other side.
    for (int r = 0; r < comm.size(); ++r) {
        for (const auto &rg : gexdss[r]) {
            auto counterpart = find_comm(gexdss[rg.dest], r);
            // Check for matching sizes
            BOOST_TEST((rg.send.size() == counterpart.recv.size()));
            BOOST_TEST((rg.recv.size() == counterpart.send.size()));

            // Check send and receive site for inconsistencies in
            // receive/send numbering.
            for (size_t i1 = 0; i1 < rg.send.size(); ++i1) {
                for (size_t i2 = i1 + 1; i2 < rg.send.size(); ++i2) {
                    auto sc1 = rg.send[i1];
                    auto sc2 = rg.send[i2];
                    auto rc1 = counterpart.recv[i1];
                    auto rc2 = counterpart.recv[i2];

                    // Recv cells could possibly be two different cells
                    // for the same send cell, if no minimal ghost layer
                    // is implemented but a "full halo".
                    // Therefore (sc1 == sc2) <=> (rc1 == rc2) might not hold.
                    // But if the send cells are different, the recv cells for
                    // sure must be different, too.
                    BOOST_TEST(if_then(sc1 != sc2, rc1 != rc2));
                }
            }

            // TODO: Cannot check contents, since indices are process local.
            //       Need ability to transform local cells to unique (global)
            //       id.
        }
    }
}

BOOST_AUTO_TEST_CASE(test_ghost_exchange_volume)
{
    boost::mpi::environment env;
    boost::mpi::communicator comm;

    for (const auto gt : repa::supported_grid_types()) {
        if (comm.rank() == 0) {
            std::cout << "Checking grid '" << repa::grid_type_to_string(gt)
                      << "'" << std::endl;
        }
        auto up = repa::make_pargrid(gt, comm, {{20., 20., 20.}}, 1.0);
        test(comm, up.get());
    }
}
