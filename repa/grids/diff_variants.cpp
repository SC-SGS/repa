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

#include "diff_variants.hpp"
#include <cassert>
#include <functional>
#include <map>
#include <memory>
#include <mpi.h>
#include <numeric>
#include <vector>

#include "util/mpi_graph.hpp"

namespace repa {
namespace grids {
namespace diff_variants {

PerNeighbor<double>
WLMVolumeComputation::compute_flow(boost::mpi::communicator neighcomm,
                                   const std::vector<rank_type> &neighbors,
                                   double load) const
{
    int nneigh = repa::util::mpi_undirected_neighbor_count(neighcomm);
    // Exchange load in local neighborhood
    std::vector<double> neighloads(nneigh);
    MPI_Neighbor_allgather(&load, 1, MPI_DOUBLE, neighloads.data(), 1,
                           MPI_DOUBLE, neighcomm);

    double avgload
        = std::accumulate(std::begin(neighloads), std::end(neighloads), load)
          / (nneigh + 1);

    // Return empty send volume if this process is underloaded
    if (load < avgload)
        return std::vector<double>(neighloads.size(), 0.0);

    std::vector<double> deficiency(neighloads.size());

    // Calculate deficiency
    for (size_t i = 0; i < neighloads.size(); ++i) {
        deficiency[i] = std::max(avgload - neighloads[i], 0.0);
    }

    auto total_deficiency
        = std::accumulate(std::begin(deficiency), std::end(deficiency), 0.0);
    double overload = load - avgload;

    // Make "deficiency" relative and then scale it to be an
    // absolute part of this process's overload
    for (size_t i = 0; i < neighloads.size(); ++i) {
        deficiency[i] = overload * deficiency[i] / total_deficiency;
    }

    return deficiency;
}

PerNeighbor<double>
SchornVolumeComputation::compute_flow(boost::mpi::communicator neighcomm,
                                      const std::vector<rank_type> &neighbors,
                                      double load) const
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

    for (uint32_t i = 0; i < _nflow_iter; i++) {
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

    return deficiency;
}

void SchornVolumeComputation::set_n_flow_iter(uint32_t nflow_iter)
{
    _nflow_iter = nflow_iter;
}

PerNeighbor<double>
SOVolumeComputation::compute_flow(boost::mpi::communicator neighcomm,
                                  const std::vector<rank_type> &neighbors,
                                  double load) const
{
    int nneigh = util::mpi_undirected_neighbor_count(neighcomm);

    std::vector<double> deficiency(nneigh);

    std::vector<int> dneigh(nneigh);
    std::vector<double> alpha(nneigh);

    MPI_Neighbor_allgather(&nneigh, 1, MPI_INT, dneigh.data(), 1, MPI_INT,
                           neighcomm);

    for (int i = 0; i < nneigh; i++)
        alpha[i] = 1.0 / (std::max(nneigh, dneigh[i]) + 1.0);

    // Exchange load in local neighborhood
    std::vector<double> neighloads(nneigh);
    MPI_Neighbor_allgather(&load, 1, MPI_DOUBLE, neighloads.data(), 1,
                           MPI_DOUBLE, neighcomm);

    if (_prev_deficiency.size() == 0) {
        for (int j = 0; j < neighloads.size(); j++)
            deficiency[j] = alpha[j] * (load - neighloads[j]);

        _prev_deficiency.reserve(nneigh);
        for (int i = 0; i < nneigh; i++)
            _prev_deficiency[neighbors[i]] = deficiency[i];
    }
    else {
        for (int j = 0; j < neighloads.size(); j++)
            deficiency[j] = _beta * alpha[j] * (load - neighloads[j])
                            + (1 - _beta) * _prev_deficiency[neighbors[j]];

        for (int i = 0; i < nneigh; i++)
            _prev_deficiency[neighbors[i]] = deficiency[i];
    }

    return deficiency;
}

void SOVolumeComputation::set_beta_value(double beta_value)
{
    _beta = beta_value;
}

PerNeighbor<double>
SOFVolumeComputation::compute_flow(boost::mpi::communicator neighcomm,
                                   const std::vector<rank_type> &neighbors,
                                   double load) const
{
    int nneigh = util::mpi_undirected_neighbor_count(neighcomm);

    std::vector<double> deficiency(nneigh);

    std::vector<int> dneigh(nneigh);
    std::vector<double> alpha(nneigh);

    MPI_Neighbor_allgather(&nneigh, 1, MPI_INT, dneigh.data(), 1, MPI_INT,
                           neighcomm);

    for (int i = 0; i < nneigh; i++)
        alpha[i] = 1.0 / (std::max(nneigh, dneigh[i]) + 1.0);

    // Exchange load in local neighborhood
    std::vector<double> neighloads(nneigh);
    MPI_Neighbor_allgather(&load, 1, MPI_DOUBLE, neighloads.data(), 1,
                           MPI_DOUBLE, neighcomm);

    std::vector<double> prev_deficiency(nneigh);
    for (int i = 0; i < _nflow_iter; i++) {
        if (i == 0) {
            for (int j = 0; j < neighloads.size(); j++)
                deficiency[j] = alpha[j] * (load - neighloads[j]);

            for (int j = 0; j < nneigh; j++)
                prev_deficiency[j] = deficiency[j];
        }
        else {
            for (int j = 0; j < neighloads.size(); j++)
                deficiency[j] = _beta * alpha[j] * (load - neighloads[j])
                                + (1 - _beta) * prev_deficiency[j];

            for (int j = 0; j < nneigh; j++)
                prev_deficiency[j] = deficiency[j];
        }
    }

    return deficiency;
}

void SOFVolumeComputation::set_n_flow_iter(uint32_t nflow_iter)
{
    _nflow_iter = nflow_iter;
}

void SOFVolumeComputation::set_beta_value(double beta_value)
{
    _beta = beta_value;
}

typedef std::function<FlowCalculator *(void)> FlowCreateFunction;
static const std::map<FlowCalcKind, FlowCreateFunction> flow_create_function_map
    = {
        {FlowCalcKind::WILLEBEEK, []() { return new WLMVolumeComputation(); }},
        {FlowCalcKind::SCHORN, []() { return new SchornVolumeComputation(); }},
        {FlowCalcKind::SO, []() { return new SOVolumeComputation(); }},
        {FlowCalcKind::SOF, []() { return new SOFVolumeComputation(); }},
};

std::unique_ptr<FlowCalculator> create_flow_calc(FlowCalcKind kind)
{
    return std::unique_ptr<FlowCalculator>(flow_create_function_map.at(kind)());
}

} // namespace diff_variants
} // namespace grids

} // namespace repa
