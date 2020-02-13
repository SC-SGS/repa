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
    int nneigh_include_center = nneigh + 1;

    // MATRIX CALCULATION
    if (_M.size() == 0) {
        _M.resize(nneigh_include_center);
        for (int i = 0; i < nneigh_include_center; i++)
            _M[i].resize(nneigh_include_center);

        for (int i = 0; i < nneigh_include_center; i++) {
            for (int j = 0; j < nneigh_include_center; j++) {
                if (i == j)
                    _M[i][j] = 0;
                else
                    _M[i][j] = (1.0 / nneigh_include_center);
            }
        }
    }

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

    // GATHERV FOR LOADS, DOES INCLUDE MYSELF
    std::vector<double> neighloads_include_center(nneigh_include_center);
    neighloads_include_center[0] = load;

    MPI_Neighbor_allgather(&load, 1, MPI_DOUBLE,
                           neighloads_include_center.data() + 1, 1, MPI_DOUBLE,
                           neighcomm);
    std::vector<int> world_load_rcounts(world_size);
    MPI_Gather(&nneigh_include_center, 1, MPI_INT, world_load_rcounts.data(), 1,
               MPI_INT, 0, comm_cart);
    int max_nneigh_include_center = *std::max_element(world_load_rcounts.begin(),
                                       world_load_rcounts.end());

    int world_load_stride = max_nneigh_include_center + 10;
    std::vector<double> world_load(world_size * world_load_stride);
    std::vector<int> world_load_displs(world_size);

    for (int i = 0; i < world_size; i++) {
        world_load_displs[i] = i * world_load_stride;
    }

    MPI_Gatherv(neighloads_include_center.data(),
                neighloads_include_center.size(), MPI_DOUBLE, world_load.data(),
                world_load_rcounts.data(), world_load_displs.data(), MPI_DOUBLE,
                0, comm_cart);
    // GATHERV FOR LOADS DONE

    // CREATE SEND_BUFFER FOR MPI_SCATTERV
    int send_buffer_stride = (*std::max_element(all_neighbors_rcounts.begin(),
                                                all_neighbors_rcounts.end()))
                             + 10;
    std::vector<double> send_buffer(world_size * send_buffer_stride);
    std::vector<int> send_buffer_counts(world_size);
    std::vector<int> send_buffer_displs(world_size);
    for (int i = 0; i < world_size; i++) {
        send_buffer_counts[i] = all_neighbors_rcounts[i];
        send_buffer_displs[i] = i * send_buffer_stride;
    }

    // CALCULATE ON RANK 0
    if (comm_cart.rank() == 0) {
        std::vector<std::vector<double>> next_load(world_size);
        if (_prev_load.size() == 0) {
            _prev_load.reserve(world_size);
            for (int j = 0; j < world_size; j++) {
                std::vector<double> w = construct_local_w(
                    world_load, world_load_rcounts, world_load_displs,
                    all_neighbors, all_neighbors_rcounts, all_neighbors_displs,
                    j);

                next_load[j] = multiply(_M, w);
            }
            _prev_load = next_load;
        }
        else {
            for (int j = 0; j < world_size; j++) {
                std::vector<double> w = construct_local_w(
                    world_load, world_load_rcounts, world_load_displs,
                    all_neighbors, all_neighbors_rcounts, all_neighbors_displs,
                    j);

                auto new_M = Matrix_scalar(_beta, _M);
                auto p1 = multiply(new_M, w);

                auto p2 = scalar((1 - _beta), _prev_load[j]);

                next_load[j] = addition(p1, p2);
            }
            _prev_load = next_load;
        }

        // MAP TO SEND_BUFFER
        for (int i = 0; i < next_load.size(); i++) {
            assert(send_buffer_counts == next_load[i]);
            int l = send_buffer_displs[i];
            for (int j = 0; j < send_buffer_counts[i]; j++) {
                send_buffer[l + j] = next_load[i][j];
            }
        }
    }

    // MPI_SCATTERV SEND_BUFFER
    std::vector<double> new_local_load(neighbors.size());
    MPI_Scatterv(send_buffer.data(), send_buffer_counts.data(),
                 send_buffer_displs.data(), MPI_DOUBLE, new_local_load.data(),
                 neighbors.size(), MPI_DOUBLE, 0, comm_cart);

    // MAP BACK TO NEIGHBORS SORTING
    std::vector<double> deficiency(neighbors.size());

    std::vector<rank_type> sorted_neighbors = neighbors;
    std::sort(sorted_neighbors.begin(), sorted_neighbors.end());
    std::unordered_map<rank_type, double> mapped_new_loads(neighbors.size());
    for (int i = 0; i < mapped_new_loads.size(); i++) {
        mapped_new_loads[sorted_neighbors[i]] = new_local_load[i];
    }

    // neighloads_include_center i + 1 to exclude ourself which is placed at
    // index 0
    for (int i = 0; i < neighbors.size(); i++) {
        // should be the other way arround but this results only in negative
        // deficiencies
        deficiency[i]
            = neighloads_include_center[i + 1] - mapped_new_loads[neighbors[i]];
    }

    return deficiency;
}

void SOCVolumeComputation::set_beta_value(double beta_value)
{
    _beta = beta_value;
}

std::vector<double> SOCVolumeComputation::construct_local_w(
    const std::vector<double> &world_load,
    const std::vector<int> &world_load_rcounts,
    const std::vector<int> &world_load_displs,
    const std::vector<rank_type> &all_neighbors,
    const std::vector<int> &all_neighbors_rcounts,
    const std::vector<int> &all_neighbors_displs,
    int j) const
{
    auto l_loads = world_load.begin() + (world_load_displs[j]);
    auto r_loads
        = world_load.begin() + (world_load_displs[j] + world_load_rcounts[j]);
    std::vector<double> w_temp(l_loads, r_loads);

    auto l_neigh = all_neighbors.begin() + (all_neighbors_displs[j]);
    auto r_neigh = all_neighbors.begin()
                   + (all_neighbors_displs[j] + all_neighbors_rcounts[j]);
    std::vector<rank_type> neighbors(l_neigh, r_neigh);

    int nneigh = neighbors.size();
    int nneigh_include_center = nneigh + 1;

    std::vector<std::tuple<rank_type, double>> w_sorted(nneigh);
    for (int i = 0; i < w_sorted.size(); i++) {
        w_sorted[i] = std::make_tuple(neighbors[i], w_temp[i + 1]);
    }
    std::sort(w_sorted.begin(), w_sorted.end());

    std::vector<double> w(nneigh_include_center);
    w[0] = w_temp[0];
    for (int i = 1; i < w.size(); i++)
        w[i] = std::get<1>(w_sorted[i - 1]);

    return w;
}

std::vector<double>
SOCVolumeComputation::addition(const std::vector<double> &v1,
                               const std::vector<double> &v2) const
{
    // get min and max size. Run array until min size and fill up array with
    // elements from vx. vx vector with more elements
    int minSize = std::min(v1.size(), v2.size());
    int maxSize = std::max(v1.size(), v2.size());

    std::vector<double> result(maxSize);
    for (int i = 0; i < minSize; i++)
        result[i] = v1[i] + v2[i];

    if (minSize != maxSize) {
        if (v1.size() > v2.size()) {
            for (int i = minSize; i < maxSize; i++)
                result[i] = v1[i];
        }
        else {
            for (int i = minSize; i < maxSize; i++)
                result[i] = v2[i];
        }
    }

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
    // use v.size() as for condition to make sure the calculations are correct.
    // Basicly shrink the matrix in x and y direction
    std::vector<double> result(v.size());

    for (int i = 0; i < v.size(); i++) {
        double r = 0.;
        for (int j = 0; j < v.size(); j++) {
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
