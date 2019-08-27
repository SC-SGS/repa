/**
 * Copyright 2017-2019 Steffen Hirschmann
 * Copyright 2017-2018 Maximilian Wildbrett
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

//#ifdef HAVE_METIS

#include "globox.hpp"
#include "glomethod.hpp"
#include "pargrid.hpp"
#include <array>
#include <parmetis.h>
#include <unordered_map>

namespace repa {
namespace grids {

struct Graph : public GloMethod {
    Graph(const boost::mpi::communicator &comm,
          Vec3d box_size,
          double min_cell_size);
    ~Graph();

private:
    friend struct HybridGPDiff; // Needs access to "partition" vector
    virtual bool sub_repartition(CellMetric m, CellCellMetric ccm) override;
};
} // namespace grids
} // namespace repa

//#endif // HAVE_METIS
