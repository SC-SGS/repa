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
    const auto c = repa::local_cell_index_type{0};
    const int which_neighbor = 0; // 0 - 26

    const auto neigh = grid->cell_neighbor_index(c, which_neighbor);
    std::cout << "Neighbor " << which_neighbor << " of cell " << c << " is ";
    if (neigh.is<repa::local_cell_index_type>())
        std::cout << " local cell " << neigh.as<repa::local_cell_index_type>();
    else
        std::cout << " ghost cell " << neigh.as<repa::ghost_cell_index_type>();
    std::cout << std::endl;
}
