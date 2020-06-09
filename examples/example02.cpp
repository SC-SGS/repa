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

/** Stream operator to output a repa::Vec3<T>.
 */
template <typename T>
std::ostream &operator<<(std::ostream &os, repa::Vec3<T> v)
{
    os << "{" << v[0] << ", " << v[1] << ", " << v[2] << "}";
    return os;
}

int main()
{
    boost::mpi::environment env;
    boost::mpi::communicator comm;

    auto grid
        = repa::make_pargrid(repa::GridType::CART, comm, {1., 1., 1.}, 1e-1);
    repa::Vec3i grid_size = grid->grid_size();
    repa::Vec3d cell_size = grid->cell_size();
    std::cout << "Grid size: " << grid_size << std::endl;
    std::cout << "Cell size: " << cell_size << std::endl;
}
