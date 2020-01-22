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

#include <memory>
#include <mpi.h>
#include <vector>

#include "pargrid.hpp"

#include "util/mpi_graph.hpp"

namespace repa {
namespace grids {
namespace diff_variants {

/**
 * 
 * Flow calculation
 * 
 */

template <typename T>
using PerNeighbor = std::vector<T>;

struct FlowCalculator {
    virtual PerNeighbor<double> compute_flow(MPI_Comm, double) const = 0;
};

struct FlowIterSetter {
    virtual void set_n_flow_iter(uint32_t nflow_iter) = 0;
};

#define DIFFUSION_MAYBE_SET_NFLOW_ITER(obj, nflow_iter)                        \
    dynamic_cast<repa::grids::diff_variants::FlowIterSetter *>(obj)            \
        ->set_n_flow_iter(nflow_iter)

/*
 * Determines the status of each process (underloaded, overloaded)
 * in the neighborhood given the local load and returns the volume of load
 * to send to each neighbor. On underloaded processes, returns a vector of
 * zeros.
 *
 * This call is collective on neighcomm.
 *
 * Default implementation follows [Willebeek Le Mair and Reeves, IEEE Tr.
 * Par. Distr. Sys. 4(9), Sep 1993] propose
 *
 * @param neighcomm Graph communicator which reflects the neighbor
 * relationship amongst processes (undirected edges), without edges to the
 *                  process itself.
 * @param load The load of the calling process.
 * @returns Vector of load values ordered according to the neighborhood
 *          ordering in neighcomm.
 */
struct WLMVolumeComputation : public FlowCalculator {
    PerNeighbor<double> compute_flow(MPI_Comm neighcomm,
                                     double load) const override;
};

struct SchornVolumeComputation : public FlowCalculator, public FlowIterSetter {
    virtual PerNeighbor<double> compute_flow(MPI_Comm neighcomm,
                                             double load) const override;
    virtual void set_n_flow_iter(uint32_t nflow_iter) override;

protected:
    uint32_t _nflow_iter = 1;
};

enum class FlowCalcKind { WILLEBEEK, SCHORN };
std::unique_ptr<FlowCalculator> create_flow_calc(FlowCalcKind);

} // namespace diff_variants
} // namespace grids

} // namespace repa
