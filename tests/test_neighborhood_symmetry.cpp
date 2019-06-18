/**
 * Checks the full shell neighborhood symmetry for inner cells.
 */

#define BOOST_TEST_MODULE neighborhood_symmetry
#include <boost/test/included/unit_test.hpp>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <cassert>
#include <repa/repa.hpp>

static bool has_neighbor(repa::grids::ParallelLCGrid *grid,
                         repa::grids::lidx d,
                         repa::grids::lidx c)
{
    for (int j = 0; j < 27; ++j) {
        if (grid->cell_neighbor_index(d, j) == c)
            return true;
    }
    return false;
}

static void test(repa::grids::ParallelLCGrid *grid)
{
    for (int c = 0; c < grid->n_local_cells(); ++c) {
        for (int j = 0; j < 27; ++j) {
            int d = grid->cell_neighbor_index(c, j);

            // If "d" is a inner cell and neighbors "c", then "c" must also
            // neighbor "d".
            if (d < grid->n_local_cells()) {
                BOOST_TEST(has_neighbor(grid, d, c));
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_neighborhood_symmetry)
{
    boost::mpi::environment env;
    boost::mpi::communicator comm;

    for (const auto gt : repa::supported_grid_types()) {
        if (comm.rank() == 0) {
            std::cout << "Checking grid '" << repa::grid_type_to_string(gt)
                      << "'" << std::endl;
        }
        auto up = repa::make_pargrid(gt, comm, {{20., 20., 20.}}, 1.0);
        test(up.get());
    }
}
