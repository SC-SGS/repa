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

#define BOOST_TEST_MODULE grid_types
#include <boost/test/included/unit_test.hpp>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include <repa/repa.hpp>

/**
 * Relative distance between a and b.
 * Only meaningfully defined for floating point types.
 */
template <typename T,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
T relative_distance(T a, T b)
{
    return std::fabs((a - b) / std::min(a, b));
}

/**
 * Returns true if the relative distance between a and b is smaller than eps.
 */
template <typename T,
          typename = std::enable_if_t<std::is_floating_point<T>::value>>
bool is_close(T a, T b, T eps = T{1e-14})
{
    return relative_distance(a, b) < eps;
}

BOOST_AUTO_TEST_CASE(test_grid_types)
{
    boost::mpi::environment env;
    boost::mpi::communicator comm;

    for (const auto gt : repa::supported_grid_types()) {
        if (comm.rank() == 0) {
            std::cout << "Checking grid '" << repa::grid_type_to_string(gt)
                      << "'" << std::endl;
        }

        const repa::Vec3d box = {{20., 20., 20.}};
        const double mings = 1.0;
        std::unique_ptr<repa::grids::ParallelLCGrid> up = nullptr;
        BOOST_CHECK_NO_THROW(up = repa::make_pargrid(gt, comm, box, mings));
        BOOST_TEST(up.get() != nullptr);

        auto grid_size = up->grid_size();
        auto cell_size = up->cell_size();
        for (size_t i = 0; i < box.size(); ++i) {
            BOOST_TEST((cell_size[i] > 0.));
            BOOST_TEST((grid_size[i] > 0));
            BOOST_TEST(grid_size[i] >= mings);
            BOOST_TEST(is_close(grid_size[i] * cell_size[i], box[i]));
        }
    }

    relative_distance(1.0, 2.0);
}
