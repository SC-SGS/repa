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
#pragma once

#include "../../pargrid.hpp"
#include "../globox.hpp"
#include "linearize.hpp"
#include <cmath>

namespace repa {
namespace util {

using namespace repa::grids; // cell index types

enum InitialPartitionType { LINEAR, CARTESIAN };

namespace impl {

template <typename AssignFunc>
void init_part_linear(const globox::GlobalBox<global_cell_index_type> &gbox,
                      const boost::mpi::communicator &comm,
                      AssignFunc &&assign_cell)
{
    auto nglobalcells = gbox.ncells();

    local_cell_index_type ncells_per_proc = static_cast<local_cell_index_type>(
        std::ceil(static_cast<double>(nglobalcells) / comm.size()));

    for (global_cell_index_type i = 0; i < nglobalcells; ++i) {
        assign_cell(i, i / ncells_per_proc);
    }
}

template <typename AssignFunc>
void init_part_cartesian(const globox::GlobalBox<global_cell_index_type> &gbox,
                         const boost::mpi::communicator &comm,
                         AssignFunc &&assign_cell)
{
    auto nglobalcells = gbox.ncells();

    int dims[3] = {0, 0, 0};
    MPI_Dims_create(comm.size(), 3, dims);

    auto cellgrid = gbox.grid_size();
    Vec3i cells_per_proc;
    for (int i = 0; i < 3; ++i) {
        cells_per_proc[i] = static_cast<Vec3i::value_type>(
            std::ceil(static_cast<double>(cellgrid[i]) / dims[i]));
    }

    for (global_cell_index_type i = 0; i < nglobalcells; ++i) {
        auto cellidx = util::unlinearize(i, gbox.grid_size());
        // Transform cellidx to 3d proc coord
        for (int i = 0; i < 3; ++i)
            cellidx[i] /= cells_per_proc[i];
        int rank;
        MPI_Cart_rank(comm, cellidx.data(), &rank);
        assign_cell(i, rank);
    }
}

} // namespace impl

struct InitPartitioning {
    InitPartitioning(const globox::GlobalBox<global_cell_index_type> &gbox,
                     const boost::mpi::communicator &comm)
        : gbox(gbox), comm(comm)
    {
    }

    template <typename AssignFunc>
    InitPartitioning &operator()(InitialPartitionType pt, AssignFunc &&f)
    {
        switch (pt) {
        case InitialPartitionType::LINEAR:
            impl::init_part_linear(gbox, comm, std::forward<AssignFunc>(f));
            break;
        case InitialPartitionType::CARTESIAN:
            impl::init_part_cartesian(gbox, comm, std::forward<AssignFunc>(f));
            break;
        }
        return *this;
    }

private:
    const globox::GlobalBox<global_cell_index_type> &gbox;
    const boost::mpi::communicator &comm;
};

} // namespace util
} // namespace repa
