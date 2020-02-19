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
                                   boost::mpi::communicator comm_cart,
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
                                      boost::mpi::communicator comm_cart,
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
SOCVolumeComputation::compute_flow(boost::mpi::communicator neighcomm,
                                   boost::mpi::communicator comm_cart,
                                   const std::vector<rank_type> &neighbors,
                                   double load) const
{
    int nneigh = util::mpi_undirected_neighbor_count(neighcomm);
    int world_size;
    MPI_Comm_size(comm_cart, &world_size);

    // GATHERV FOR NEIGHBORS, DOES NOT INCLUDE MYSELF
    std::vector<int> all_neighbors_rcounts(world_size);
    MPI_Gather(&nneigh, 1, MPI_INT, all_neighbors_rcounts.data(), 1, MPI_INT, 0,
               comm_cart);
    int max_nneigh = *std::max_element(all_neighbors_rcounts.begin(),
                                       all_neighbors_rcounts.end());

    int all_neighbors_stride = max_nneigh + 10;
    std::vector<int> all_neighbors(world_size * all_neighbors_stride);
    std::vector<int> all_neighbors_displs(world_size);

    for (int i = 0; i < world_size; i++) {
        all_neighbors_displs[i] = i * all_neighbors_stride;
    }

    MPI_Gatherv(neighbors.data(), neighbors.size(), MPI_INT,
                all_neighbors.data(), all_neighbors_rcounts.data(),
                all_neighbors_displs.data(), MPI_INT, 0, comm_cart);
    // GATHERV FOR NEIGHBORS DONE

    // GATHER LOAD
    std::vector<double> w(world_size);
    MPI_Gather(&load, 1, MPI_DOUBLE, w.data(), 1, MPI_DOUBLE, 0, comm_cart);

    // CALCULATE ON RANK 0
    std::vector<double> next_load(world_size);
    if (comm_cart.rank() == 0) {
        // MATRIX CALCULATION
        if (_M.size() == 0) {
            _M.resize(world_size);
            for (int i = 0; i < _M.size(); i++)
                _M[i].resize(world_size);

            for (int i = 0; i < _M.size(); i++) {
                auto l_neigh
                    = all_neighbors.begin() + (all_neighbors_displs[i]);
                auto r_neigh
                    = all_neighbors.begin()
                      + (all_neighbors_displs[i] + all_neighbors_rcounts[i]);
                for (int j = 0; j < _M.size(); j++) {
                    if (i == j)
                        _M[i][j] = 0;
                    else if (std::find(l_neigh, r_neigh, j) != r_neigh)
                        _M[i][j] = (1.0 / all_neighbors_rcounts[i]);
                    else
                        _M[i][j] = 0;
                }
            }
        }

        if (_prev_load.size() == 0) {
            next_load = multiply(_M, w);
            _prev_load = w;
        }
        else {
            auto new_M = Matrix_scalar(_beta, _M);
            auto p1 = multiply(new_M, w);
            auto p2 = scalar((1 - _beta), _prev_load);

            next_load = addition(p1, p2);
            _prev_load = w;
        }
    }

    // MPI_Bcast next_load
    MPI_Bcast(next_load.data(), next_load.size(), MPI_DOUBLE, 0, comm_cart);

    std::vector<double> deficiency(nneigh);
    double d = load - next_load[comm_cart.rank()];
    for (int i = 0; i < deficiency.size(); i++)
        deficiency[i] = (1. / nneigh) * d;

    return deficiency;
}

void SOCVolumeComputation::set_beta_value(double beta_value)
{
    _beta = beta_value;
}

std::vector<double>
SOCVolumeComputation::addition(const std::vector<double> &v1,
                               const std::vector<double> &v2) const
{
    assert(v1.size() == v2.size());
    std::vector<double> result(v1.size());
    for (int i = 0; i < v1.size(); i++)
        result[i] = v1[i] + v2[i];

    return result;
}

std::vector<double>
SOCVolumeComputation::scalar(double scalar, const std::vector<double> &v) const
{
    std::vector<double> result(v.size());
    for (int i = 0; i < v.size(); i++)
        result[i] = scalar * v[i];

    return result;
}

std::vector<std::vector<double>> SOCVolumeComputation::Matrix_scalar(
    double scalar, const std::vector<std::vector<double>> &M) const
{
    std::vector<std::vector<double>> result(M.size(),
                                            std::vector<double>(M.size()));
    for (int i = 0; i < M.size(); i++)
        for (int j = 0; j < M[i].size(); j++)
            result[i][j] = scalar * M[i][j];

    return result;
}

std::vector<double>
SOCVolumeComputation::multiply(const std::vector<std::vector<double>> &M,
                               const std::vector<double> &v) const
{
    std::vector<double> result(v.size());

    for (int i = 0; i < M.size(); i++) {
        double r = 0.;
        for (int j = 0; j < M[i].size(); j++) {
            r += M[i][j] * v[j];
        }
        result[i] = r;
    }

    return result;
}

PerNeighbor<double>
SOVolumeComputation::compute_flow(boost::mpi::communicator neighcomm,
                                  boost::mpi::communicator comm_cart,
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
                                   boost::mpi::communicator comm_cart,
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
        {FlowCalcKind::SOC, []() { return new SOCVolumeComputation(); }},
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
