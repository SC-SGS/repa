/**
 * Copyright 2017-2019 Steffen Hirschmann
 * Copyright 2017-2018 Maximilian Wildbrett
 * Copyright 2019      Simon Hauser
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

struct Schornbaum : public Diffusion {
    Schornbaum(const boost::mpi::communicator &comm,
               Vec3d box_size,
               double min_cell_size);
    ~Schornbaum();

    void command(std::string s) override;

protected:
    /*
     * See Schornbaum, F.; Rüde, U. Extreme-Scale Block-Structured Adaptive Mesh
     * Refinement. SIAM J. Sci. Comput. 2018, 40, C358–C387.
     */
    virtual std::vector<double> compute_send_volume(double load) override;

private:
    int flow_count = 15;
};
} // namespace grids
} // namespace repa
