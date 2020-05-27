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

#include "glomethod.hpp"

namespace repa {
namespace grids {

/** Regular Cartesian process grid; equally sized boxes.
 *
 * Cells are ordered on the ghost grid according to a simple row-wise ordering
 * (see impl::linearize). All 3d indices live on the ghost grid, i.e.
 * their values range from {0, 0, 0} to m_ghost_grid.
 * Where {0, 0, 0} is the first ghost cell, {1, 1, 1} the first inner (i.e.
 * boundary) cell and so on.
 */
struct CartGrid : public GloMethod {
    CartGrid(const boost::mpi::communicator &comm,
             Vec3d box_size,
             double min_cell_size,
             ExtraParams ep);

    virtual ~CartGrid();

    util::ioptional<rank_type>
    rank_of_cell(global_cell_index_type idx) const override;

    bool sub_repartition(CellMetric m, CellCellMetric ccm) override;
};
} // namespace grids
} // namespace repa
