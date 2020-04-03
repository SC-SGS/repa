/**
 * Copyright 2017-2020 Steffen Hirschmann
 * Copyright 2020 Benjamin Vier
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
 * Checks mpi_cart_coloring.hpp
 */

#define BOOST_TEST_MODULE coloring

#include <vector>

#include <boost/mpi.hpp>
#include <boost/test/included/unit_test.hpp>
#include <repa/grids/util/mpi_cart_coloring.hpp>
#include <repa/grids/util/vec_arith.hpp>

static int neighbor_rank(const boost::mpi::communicator &comm, repa::Vec3i off)
{
    auto coord = repa::util::mpi_cart_get_coords(comm);
    const auto dims = repa::util::mpi_cart_get_dims(comm);

    using namespace repa::util::vector_arithmetic;
    coord = (coord + off) % dims;
    return repa::util::mpi_cart_rank(comm, coord);
}

static bool iff(bool a, bool b)
{
    return (a && b) || (!a && !b);
}

/**
 * Test the validity of the coloring.
 * (Check if sets are independent.)
 */
BOOST_AUTO_TEST_CASE(test_coloring)
{
    boost::mpi::environment env;
    boost::mpi::communicator comm;
    repa::Vec3i dims{0, 0, 0}, periods{1, 1, 1};
    MPI_Dims_create(comm.size(), 3, dims.data());

    MPI_Comm _comm_cart;
    MPI_Cart_create(comm, 3, dims.data(), periods.data(), true, &_comm_cart);

    boost::mpi::communicator comm_cart(_comm_cart,
                                       boost::mpi::comm_take_ownership);

    int my_color = -1;
    {
        int cur_color = 0;
        repa::util::independent_process_sets(comm_cart)
            .for_each([&]() {
                BOOST_CHECK(my_color == -1);
                my_color = cur_color;
            })
            .for_all_after_each_round([&]() { cur_color++; })();
    }
    BOOST_CHECK(my_color != -1);

    std::vector<int> colors;
    boost::mpi::all_gather(comm, my_color, colors);

    if (comm_cart.rank() == 0) {
        std::cout << "Colors:";
        for (auto c : colors)
            std::cout << " " << c;
        std::cout << std::endl;
    }

    for (int rank = 0; rank < comm.size(); ++rank) {
        repa::Vec3i off;
        for (off[0] = -1; off[0] <= 1; ++off[0]) {
            for (off[1] = -1; off[1] <= 1; ++off[1]) {
                for (off[2] = -1; off[2] <= 1; ++off[2]) {
                    auto neighrank = neighbor_rank(comm_cart, off);
                    BOOST_CHECK(iff(rank != neighrank,
                                    colors.at(rank) != colors.at(neighrank)));
                }
            }
        }
    }
}
