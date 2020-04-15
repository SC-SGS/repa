/**
 * Copyright 2017-2020 The repa authors
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

#include <boost/mpi/communicator.hpp>
#include <memory>
#include <mpi.h>
#include <unordered_map>
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
    /*
     * Determines the status of each process (underloaded, overloaded)
     * in the neighborhood given the local load and returns the volume of load
     * to send to each neighbor. On underloaded processes, returns a vector of
     * zeros.
     *
     * This call is collective on neighcomm.
     *
     * @param neighcomm Graph communicator which reflects the neighbor
     * relationship amongst processes (undirected edges), without edges to the
     *                  process itself.
     * @param neighbors The ranks of the neighbors of the calling process.
     * @param load The load of the calling process.
     * @returns Vector of load values ordered according to the neighborhood
     *          ordering in neighcomm.
     */
    virtual PerNeighbor<double>
    compute_flow(boost::mpi::communicator neighcomm,
                 const std::vector<rank_type> &neighbors,
                 double load) const = 0;
};

struct FlowIterSetter {
    virtual void set_n_flow_iter(uint32_t nflow_iter) = 0;
};

struct BetaValueSetter {
    virtual void set_beta_value(double beta_value) = 0;
};

inline bool diffusion_maybe_set_nflow_iter(FlowCalculator *ptr,
                                           uint32_t nflow_iter)
{
    FlowIterSetter *f;
    if ((f = dynamic_cast<FlowIterSetter *>(ptr)))
        f->set_n_flow_iter(nflow_iter);
    return f != nullptr;
}

inline bool diffusion_maybe_set_beta(FlowCalculator *ptr, double beta)
{
    BetaValueSetter *b;
    if ((b = dynamic_cast<BetaValueSetter *>(ptr)))
        b->set_beta_value(beta);
    return b != nullptr;
}

/*
 * This implementation follows [Willebeek Le Mair and Reeves, IEEE Tr. Par.
 * Distr. Sys. 4(9), Sep 1993] propose.
 */
struct WLMVolumeComputation : public FlowCalculator {
    virtual PerNeighbor<double>
    compute_flow(boost::mpi::communicator neighcomm,
                 const std::vector<rank_type> &neighbors,
                 double load) const override;
};

/*
 * This implementation follows [Florian Schornbaum and Ulrich Rüde, SIAM J. Sci.
 * Comput., 40(3), C358–C387.] propose.
 *
 * Flow will be set with set_n_flow_iter(uint32_t nflow_iter), by using
 * dd.command("set flow_count 15")
 * Flow Default is 1. (results in using [Cybenko, Journal of Parallel and
 * Distributed Computing Volume 7, Issue 2, October 1989, Pages 279-301]
 * propose.)
 */
struct SchornVolumeComputation : public FlowCalculator, public FlowIterSetter {
    virtual PerNeighbor<double>
    compute_flow(boost::mpi::communicator neighcomm,
                 const std::vector<rank_type> &neighbors,
                 double load) const override;
    virtual void set_n_flow_iter(uint32_t nflow_iter) override;

protected:
    uint32_t _nflow_iter = 1;
};

/*
 * This implementation follows [Muthukrishnan, Ghosh and Schultz, Theory of
 * Computing Systems volume 31, pages331–354(1998)] propose.
 *
 * Beta will be set with set_beta_value(double beta_value), by using
 * dd.command("set beta <value>")
 * Beta Default: 1.8
 */
struct SOVolumeComputation : public FlowCalculator, public BetaValueSetter {
    virtual PerNeighbor<double>
    compute_flow(boost::mpi::communicator neighcomm,
                 const std::vector<rank_type> &neighbors,
                 double load) const override;

    virtual void set_beta_value(double beta_value) override;

protected:
    double _beta = 1.8;

private:
    mutable std::unordered_map<rank_type, double> _prev_deficiency;
};

/*
 * This implementation follows [Muthukrishnan, Ghosh and Schultz, Theory of
 * Computing Systems volume 31, pages331–354(1998)] propose with a minor change.
 *
 * Beta will be set with set_beta_value(double beta_value), by using
 * dd.command("set beta <value>")
 * Beta Default: 1.8
 *
 * The change is that in each repartition step the flow calculation can be
 * iterated by setting the flow variable using the method set_n_flow_iter.
 * This method will be called when using dd.command("set flow_count 15").
 * Flow Default is 1.
 */
struct SOFVolumeComputation : public FlowCalculator,
                              public FlowIterSetter,
                              public BetaValueSetter {
    virtual PerNeighbor<double>
    compute_flow(boost::mpi::communicator neighcomm,
                 const std::vector<rank_type> &neighbors,
                 double load) const override;

    virtual void set_n_flow_iter(uint32_t nflow_iter) override;
    virtual void set_beta_value(double beta_value) override;

protected:
    double _beta = 1.8;
    uint32_t _nflow_iter = 1;

private:
    mutable std::unordered_map<rank_type, double> _prev_deficiency;
};

enum class FlowCalcKind { WILLEBEEK, SCHORN, SO, SOF };
std::unique_ptr<FlowCalculator> create_flow_calc(FlowCalcKind);

} // namespace diff_variants
} // namespace grids

} // namespace repa
