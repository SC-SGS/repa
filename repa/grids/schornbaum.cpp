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

#include "schornbaum.hpp"
#include <regex>

#include "util/mpi_graph.hpp"

#ifndef NDEBUG
#define SCHORNBAUM_DEBUG
#endif

namespace repa {
namespace grids {

std::vector<double> Schornbaum::compute_send_volume(double load)
{
    int nneigh = util::mpi_undirected_neighbor_count(neighcomm);
    // Init flow(deficiency) and alpha
    std::vector<double> deficiency(nneigh);

    std::vector<int> dneigh(nneigh);
    std::vector<double> alpha(nneigh);

    MPI_Neighbor_allgather(&nneigh, 1, MPI_INT, dneigh.data(), 1, MPI_INT,
                           neighcomm);

    for (int i = 0; i < nneigh; i++)
        alpha[i] = 1.0 / (std::max(nneigh, dneigh[i]) + 1.0);

#ifdef SCHORNBAUM_DEBUG
    int world_size;
    MPI_Comm_size(comm, &world_size);
    std::vector<double> globalLoadListBefore(world_size);
    MPI_Allgather(&load, 1, MPI_DOUBLE, globalLoadListBefore.data(), 1,
                  MPI_DOUBLE, comm);
    double globalLoadSumBefore = std::accumulate(
        std::begin(globalLoadListBefore), std::end(globalLoadListBefore), 0.0);
#endif

    for (int i = 0; i < flow_count; i++) {
        // Exchange load in local neighborhood
        std::vector<double> neighloads(nneigh);
        MPI_Neighbor_allgather(&load, 1, MPI_DOUBLE, neighloads.data(), 1,
                               MPI_DOUBLE, neighcomm);

        double old_load = load;
        for (int j = 0; j < neighloads.size(); j++) {
            double new_f = alpha[j] * (old_load - neighloads[j]);
            deficiency[j] += new_f;
            load -= new_f;
        }
    }

#ifdef SCHORNBAUM_DEBUG
    std::vector<double> globalLoadListAfter(world_size);
    MPI_Allgather(&load, 1, MPI_DOUBLE, globalLoadListAfter.data(), 1,
                  MPI_DOUBLE, comm);

    double globalLoadSumAfter = std::accumulate(
        std::begin(globalLoadListAfter), std::end(globalLoadListAfter), 0.0);

    assert(std::abs(globalLoadSumBefore - globalLoadSumAfter) <= 0.001);
#endif

    return deficiency;
}

/*
 * Initialization
 */
Schornbaum::Schornbaum(const boost::mpi::communicator &comm,
                       Vec3d box_size,
                       double min_cell_size)
    : Diffusion(comm, box_size, min_cell_size)
{
}

Schornbaum::~Schornbaum()
{
}

void Schornbaum::command(std::string s)
{
    // Example: set flow_count 15
    static const std::regex mure("(set) (flow_count) ([[:digit:]]+)");
    std::smatch m;

    if (std::regex_match(s, m, mure)) {
        flow_count = std::stoi(m[3].str().c_str(), NULL);
        std::cout << "Setting flow_count = " << flow_count << std::endl;
    }
}
} // namespace grids
} // namespace repa
