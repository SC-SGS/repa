/**
 * Copyright 2017-2020 Steffen Hirschmann
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

#include <boost/mpi.hpp>
#include <iostream>
#include <repa/repa.hpp>

int main()
{
    boost::mpi::environment env;
    boost::mpi::communicator comm;

    auto grid
        = repa::make_pargrid(repa::GridType::CART, comm, {1., 1., 1.}, 1e-1);

    // Use the range objects in conjunction with range-based for loops
    for (const auto local_idx : grid->local_cells()) {
        // Index type cast implicitly to int.
        std::cout << "Rank " << comm.rank()
                  << " local cell index: " << local_idx;
        // On debug builds the function "global_hash" will return the associated
        // global cell index to a local index:
        std::cout << " global index: " << grid->global_hash(local_idx)
                  << std::endl;
    }
}
