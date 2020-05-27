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

#include "cart.hpp"
#include "util/linearize.hpp"
#include "util/mpi_cart.hpp"
#include "util/vec_arith.hpp"

namespace repa {
namespace grids {

/** FIXME: Build atop box_global_cell_index_storage
 */

CartGrid::CartGrid(const boost::mpi::communicator &comm,
                   Vec3d box_size,
                   double min_cell_size,
                   ExtraParams ep)
    : GloMethod(comm, box_size, min_cell_size, ep)
{
}

CartGrid::~CartGrid()
{
}

util::ioptional<rank_type>
CartGrid::rank_of_cell(global_cell_index_type idx) const
{
    using namespace util::vector_arithmetic;
    const auto idx3d = util::unlinearize(idx, gbox.grid_size());
    const auto dims = util::mpi_cart_get_dims(comm_cart);
    const Vec3d cells_per_proc
        = static_cast_vec<Vec3d>(gbox.grid_size()) / dims;

    Vec3d tmp = idx3d / cells_per_proc - .5;
    for (size_t i = 0; i < tmp.size(); ++i)
        tmp[i] = std::round(tmp[i]);

    Vec3i procidx
        = vec_clamp(static_cast_vec<Vec3i>(tmp), constant_vec3(0), dims - 1);
    return util::mpi_cart_rank(comm_cart, procidx);
}

bool CartGrid::sub_repartition(CellMetric m, CellCellMetric ccm)
{
    UNUSED(m);
    UNUSED(ccm);
    return false;
}

} // namespace grids
} // namespace repa
