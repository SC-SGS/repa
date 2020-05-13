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

#include <algorithm>
#include <cctype>
#include <vector>

#include "initial_partitioning.hpp"

namespace repa {
namespace util {

namespace {

static const struct {
    std::vector<std::string> descriptors;
    InitialPartitionType pt;
} partition_registry[] = {
    {{"linear"}, InitialPartitionType::LINEAR},
    {{"cart1d", "cartesian1d"}, InitialPartitionType::CARTESIAN1D},
    {{"cart3d", "cartesian3d", "cart"}, InitialPartitionType::CARTESIAN3D},
};

} // namespace

util::InitialPartitionType parse_part_type(const std::string &desc)
{
    // Ignore case
    std::string ldesc;
    ldesc.reserve(desc.size());
    std::transform(std::begin(desc), std::end(desc), std::back_inserter(ldesc),
                   [](char c) { return std::tolower(c); });

    for (const auto &el : partition_registry) {
        if (std::find(std::begin(el.descriptors), std::end(el.descriptors),
                      ldesc)
            != std::end(el.descriptors)) {
            return el.pt;
        }
    }

    throw UnknownInitialPartitionTypeError(desc);
}

boost::mpi::communicator
make_init_part_communicator(const boost::mpi::communicator &comm,
                            InitialPartitionType init_part)
{
    switch (init_part) {
    case util::InitialPartitionType::CARTESIAN1D:
        return boost::mpi::communicator{
            make_cartesian_communicator(comm, Vec3i{0, 1, 1}),
            boost::mpi::comm_take_ownership};
        break;
    case util::InitialPartitionType::CARTESIAN3D:
        return boost::mpi::communicator{make_cartesian_communicator(comm),
                                        boost::mpi::comm_take_ownership};
        break;
    default:
        return boost::mpi::communicator{MPI_COMM_NULL, boost::mpi::comm_attach};
        break;
    }
}

} // namespace util
} // namespace repa
