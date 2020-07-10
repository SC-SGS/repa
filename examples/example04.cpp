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
    const auto pos = repa::Vec3d{.75, .75, .75};

    const int rank = grid->position_to_rank(pos);
    if (rank == comm.rank()) {
        // On ranks that do not own "pos" the following call will throw
        const auto cellidx = grid->position_to_cell_index(pos);
        std::cout
            << pos << " belongs to rank "
            << rank
            // "cellidx" is of type repa::local_cell_index_type and can be
            // (implicitly) converted to an integer representing the cell index.
            << " into cell: " << static_cast<int>(cellidx) << std::endl;
    }
}
