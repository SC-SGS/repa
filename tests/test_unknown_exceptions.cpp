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
#define BOOST_TEST_MODULE unknown_exceptions
#include <boost/test/unit_test.hpp>

#include <boost/mpi.hpp>
#include <repa/repa.hpp>

template <typename E>
bool ignore_message(const E &e)
{
    (void)e;
    return true;
}

bool message_unknown_grid_type(const repa::UnknownGridTypeError &e)
{
    return std::string{e.what()} == "Unknown grid type: unknown";
}

bool message_unknown_grid_type_ff(const repa::UnknownGridTypeError &e)
{
    return std::string{e.what()} == "Unknown grid type: 255";
}

bool message_unknown_command(
    const repa::grids::ParallelLCGrid::UnknwonCommandError &e)
{
    return std::string{e.what()} == "Unknown command: unknown";
}

BOOST_AUTO_TEST_CASE(test_unknown_exceptions)
{
    repa::GridType gt;
    BOOST_CHECK_EXCEPTION(gt = repa::parse_grid_type("unknown"),
                          repa::UnknownGridTypeError,
                          message_unknown_grid_type);

    // Some invalid grid type
    gt = static_cast<repa::GridType>(0xFF);
    BOOST_CHECK_EXCEPTION(repa::has_grid_type(gt), repa::UnknownGridTypeError,
                          ignore_message);
    BOOST_CHECK_EXCEPTION(repa::grid_type_to_string(gt),
                          repa::UnknownGridTypeError, ignore_message);

    boost::mpi::communicator comm;
    BOOST_CHECK_EXCEPTION(repa::make_pargrid(gt, comm, {1., 1., 1.}, 1.),
                          repa::UnknownGridTypeError,
                          message_unknown_grid_type_ff);

    gt = repa::GridType::CART;
    auto g = repa::make_pargrid(gt, comm, {1., 1., 1.}, 1. / (comm.size()));
    BOOST_CHECK_EXCEPTION(g->command("unknown"),
                          repa::grids::ParallelLCGrid::UnknwonCommandError,
                          message_unknown_command);
}

int main(int argc, char **argv)
{
    boost::mpi::environment mpi_env{argc, argv};
    return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
