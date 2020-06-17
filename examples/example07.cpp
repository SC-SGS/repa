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

// Let's see how often we come by global cell index "c_want"
const int c_want = 120;
void f(repa::global_cell_index_type c, repa::global_cell_index_type d)
{
    static int ncalls = 0;
    if (c == c_want || d == c_want) {
        ncalls++;
        std::cout << ncalls << " ";
        std::cout << c << ", " << d << std::endl;
    }
}

int main()
{
    boost::mpi::environment env;
    boost::mpi::communicator comm;

    auto grid
        = repa::make_pargrid(repa::GridType::CART, comm, {1., 1., 1.}, 1e-1);
    for (const auto &c : grid->local_cells()) {
        for (int which_neigh = 0; which_neigh <= 13; ++which_neigh) {
            const auto d = grid->cell_neighbor_index(c, which_neigh);

            // Gets called exactly once for all unordered pairs {c, d} (which
            // includes the pair {c, c}). We use global cells here for
            // demonstration purposes.
            f(grid->global_hash(c), grid->global_hash(d));
        }
    }
}
