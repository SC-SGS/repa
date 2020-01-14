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

enum InitialPartitionType { LINEAR, CARTESIAN1D, CARTESIAN3D };

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

namespace __cart_impl {

constexpr inline repa::Vec3i ceil_div(const repa::Vec3i &a,
                                      const repa::Vec3i &b)
{
    Vec3i result;
    for (int i = 0; i < 3; ++i) {
        // Does not assign an equal amount of cells to each process.
        result[i] = static_cast<Vec3i::value_type>(
            std::ceil(static_cast<double>(a[i]) / b[i]));
    }
    return result;
}

struct Cart_CellProcessIndexConverter {
    constexpr Cart_CellProcessIndexConverter(const repa::Vec3i &cell_grid,
                                             const repa::Vec3i &dims)
        : dims(dims), cells_per_proc(ceil_div(cell_grid, dims))
    {
    }

    constexpr Vec3i operator()(const Vec3i &cell_idx) const
    {
        Vec3i result;
        for (int d = 0; d < 3; ++d) {
            result[d] = std::min(cell_idx[d] / cells_per_proc[d], dims[d] - 1);
        }
        return result;
    }

private:
    const repa::Vec3i &dims;
    const repa::Vec3i cells_per_proc;
};

} // namespace __cart_impl

template <typename AssignFunc>
void init_part_cartesian3d(
    const globox::GlobalBox<global_cell_index_type> &gbox,
    const boost::mpi::communicator &comm,
    AssignFunc &&assign_cell)
{
    assert(comm.has_cartesian_topology());
    const auto nglobalcells = gbox.ncells();
    const auto cellgrid = gbox.grid_size();

    Vec3i dims{0, 0, 0};
    MPI_Dims_create(comm.size(), 3, dims.data());
    const __cart_impl::Cart_CellProcessIndexConverter to_procidx{cellgrid,
                                                                 dims};

    for (global_cell_index_type i = 0; i < nglobalcells; ++i) {
        const auto procidx = to_procidx(util::unlinearize(i, cellgrid));
        int rank;
        MPI_Cart_rank(comm, procidx.data(), &rank);
        assign_cell(i, rank);
    }
}

template <typename AssignFunc>
void init_part_cartesian1d(
    const globox::GlobalBox<global_cell_index_type> &gbox,
    const boost::mpi::communicator &comm,
    AssignFunc &&assign_cell)
{
    assert(comm.has_cartesian_topology());
    const auto nglobalcells = gbox.ncells();
    const auto cellgrid = gbox.grid_size();

    Vec3i dims{comm.size(), 1, 1};
    const __cart_impl::Cart_CellProcessIndexConverter to_procidx{cellgrid,
                                                                 dims};

    for (global_cell_index_type i = 0; i < nglobalcells; ++i) {
        const auto procidx = to_procidx(util::unlinearize(i, gbox.grid_size()));
        // 1d partitioning, thus, rank is procidx[0].
        // Do not use MPI_Cart_rank here, as it is a 3d Cartesian communicator.
        assign_cell(i, procidx[0]);
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
        case InitialPartitionType::CARTESIAN1D:
            impl::init_part_cartesian1d(gbox, comm,
                                        std::forward<AssignFunc>(f));
            break;
        case InitialPartitionType::CARTESIAN3D:
            impl::init_part_cartesian3d(gbox, comm,
                                        std::forward<AssignFunc>(f));
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
