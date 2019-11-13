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

#include "repa/repa.hpp"
#include <boost/mpi/communicator.hpp>
#include <functional>
#include <set>

namespace testenv {

using MetricFunc = std::function<std::vector<double>(size_t)>;

struct TEnv {
    boost::mpi::communicator comm;
    repa::Vec3d box;
    double mings;
    repa::ExtraParams ep;
    std::set<repa::GridType> grids;
    bool repart = true, repart_twice = false;
    MetricFunc *get_metric = nullptr;

    TEnv(const boost::mpi::communicator &comm,
         repa::Vec3d box,
         double mings,
         repa::ExtraParams ep);
    TEnv(repa::ExtraParams ep);

    TEnv &with_repart();
    TEnv &with_repart(MetricFunc &f);
    TEnv &with_repart_twice();
    TEnv &without_repart();
    TEnv &all_grids();
    TEnv &only(std::set<repa::GridType> s);
    TEnv &exclude(std::set<repa::GridType> s);

    using Fno_gt
        = std::function<void(const TEnv &, repa::grids::ParallelLCGrid *)>;
    using Fwith_gt = std::function<void(
        const TEnv &, repa::grids::ParallelLCGrid *, repa::GridType)>;

   void run(Fno_gt test_func);
   void run(Fwith_gt test_func);

private:
    using TestFunc = std::function<void(
        const TEnv &te, repa::grids::ParallelLCGrid *grid, repa::GridType gt)>;

    void run_impl(TestFunc test_func);
};
} // namespace testenv

testenv::TEnv default_test_env(repa::ExtraParams ep = repa::ExtraParams{});
