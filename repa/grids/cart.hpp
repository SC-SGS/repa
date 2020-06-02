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

/** Regular Cartesian process grid; equally sized boxes, if divisible. Else
 * some processes get larger boxes than others.
 *
 * Supports all initial partitioning initializations: linear or Cartesian 1/2/3d
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

private:
    util::StaticRankAssigner<decltype(gbox)> _static_part;
};
} // namespace grids
} // namespace repa
