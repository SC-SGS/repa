/**
 * Checks for the correct creation of supported pargrids.
 */

#define BOOST_TEST_MODULE grid_types
#include <boost/test/included/unit_test.hpp>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <cassert>
#include <repa/repa.hpp>

BOOST_AUTO_TEST_CASE(test_grid_types)
{
    boost::mpi::environment env;
    boost::mpi::communicator comm;

    for (const auto gt : repa::supported_grid_types()) {
        std::cout << "Checking grid '" << repa::grid_type_to_string(gt) << "'"
                  << std::endl;
        auto up = repa::make_pargrid(gt, comm, {{10., 10., 10.}}, 1.0);
        BOOST_TEST(up.get() != nullptr);
    }
}
