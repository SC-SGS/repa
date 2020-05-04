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
#include "vec_arith.hpp"
#include <cmath>

namespace repa {
namespace util {

using namespace repa::grids; // cell index types

enum InitialPartitionType : int { LINEAR, CARTESIAN1D, CARTESIAN3D };

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

struct Cart_CellProcessIndexConverter {
    static repa::Vec3d _vdiv(const Vec3i &a, const Vec3i &b)
    {
        using namespace repa::util::vector_arithmetic;
        return static_cast_vec<Vec3d>(a) / static_cast_vec<Vec3d>(b);
    }

    static Vec3i _round_vec(const Vec3d &a)
    {
        Vec3i res;
        for (int d = 0; d < 3; ++d)
            res[d] = static_cast<int>(std::round(a[d]));
        return res;
    }

    Cart_CellProcessIndexConverter(const repa::Vec3i &cell_grid,
                                   const repa::Vec3i &dims)
        : dims(dims), cells_per_proc(_vdiv(cell_grid, dims))
    {
    }

    Vec3i operator()(const Vec3i &cell_idx) const
    {
        using namespace repa::util::vector_arithmetic;
        return vec_clamp(_round_vec(cell_idx / cells_per_proc - .5),
                         constant_vec3(0), dims - 1);
    }

private:
    const repa::Vec3i &dims;
    const repa::Vec3d cells_per_proc;
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

    const double layers_per_proc = cellgrid[0] / comm.size();

    for (global_cell_index_type i = 0; i < nglobalcells; ++i) {
        const int xindex = util::unlinearize(i, gbox.grid_size())[0];
        // Use std::round to get a more equal distribution.
        const int procidx
            = std::min(static_cast<int>(std::round(xindex / layers_per_proc)),
                       comm.size() - 1);
        assign_cell(i, procidx);
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
