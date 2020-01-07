/**
 * Copyright 2017-2019 The repa authors
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

#include "diffusion.hpp"

namespace repa {
namespace grids {

struct PSDiffusion : public Diffusion {
    PSDiffusion(const boost::mpi::communicator &comm,
                Vec3d box_size,
                double min_cell_size);
    ~PSDiffusion();

protected:
    virtual std::vector<double> compute_send_volume(double load) override;

    virtual void post_init(bool firstcall) override;

    virtual PerNeighbor<GlobalCellIndices>
    compute_send_list(std::vector<double> &&sendLoads,
                      const std::vector<double> &weights) override;

private:
    bool coords_based_allow_sending(local_cell_index_type c,
                                    rank_type neighrank);
    Vec3i fix_periodic_edge(const Vec3i &c0, const Vec3i &c2);

#ifndef NDEBUG
    bool rank_based_allow_sending(local_cell_index_type c, rank_type neighrank);

    std::vector<rank_type> neighborhood_ranks;
    std::unordered_map<int, rank_type *> nr_mappings;

    std::vector<rank_type> initial_neighborhood;
#endif
};

} // namespace grids
} // namespace repa
