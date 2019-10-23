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
#include <boost/test/included/unit_test.hpp>

#include "repa/repa.hpp"
#include <boost/mpi/communicator.hpp>
#include <functional>
#include <random>
#include <set>

inline void repartition_randomly(repa::grids::ParallelLCGrid *grid)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1., 10.);

    auto random_dbl = [&gen, &dis]() { return dis(gen); };

    auto nlc = grid->n_local_cells();
    auto random_cellweights = [&random_dbl, nlc]() {
        std::vector<double> w(nlc);
        std::generate(std::begin(w), std::end(w), random_dbl);
        return w;
    };
    auto noop = []() {};
    auto ccm = [](int, int) { return 1.0; };

    BOOST_CHECK_NO_THROW(grid->repartition(random_cellweights, ccm, noop));
}

namespace {
struct TEnv {
    boost::mpi::communicator comm;
    repa::Vec3d box;
    double mings;
    repa::ExtraParams ep;
    std::set<repa::GridType> grids;
    bool repart = true, repart_twice = false;

    TEnv(const boost::mpi::communicator &comm,
         repa::Vec3d box,
         double mings,
         repa::ExtraParams ep = repa::ExtraParams{})
        : comm(comm, boost::mpi::comm_duplicate),
          box(std::move(box)),
          mings(mings),
          ep(ep)
    {
    }

    TEnv(repa::ExtraParams ep = repa::ExtraParams{}) : mings(1.0), ep(ep)
    {
        // Devise some appropriately sized grid suitable for all methods.
        repa::Vec3i dims{0, 0, 0};
        MPI_Dims_create(comm.size(), 3, dims.data());
        for (size_t i = 0; i < box.size(); ++i)
            box[i] = mings * 5 * (dims[i] + 1);
    }

    TEnv &with_repart()
    {
        repart = true;
        repart_twice = false;
        return *this;
    }

    TEnv &with_repart_twice()
    {
        repart = true;
        repart_twice = true;
        return *this;
    }

    TEnv &without_repart()
    {
        repart = false;
        repart_twice = false;
        return *this;
    }

    TEnv &all_grids()
    {
        grids = repa::supported_grid_types();
        return *this;
    }

    TEnv &only(std::set<repa::GridType> s)
    {
        grids.clear();
        auto supported = repa::supported_grid_types();
        for (const auto gt : s) {
            if (supported.find(gt) != std::end(supported))
                grids.insert(gt);
        }
        return *this;
    }

    TEnv &exclude(std::set<repa::GridType> s)
    {
        for (const auto gt : s)
            grids.erase(gt);
        return *this;
    }

    template <typename Func>
    void run(Func test_func)
    {
        for (const auto gt : grids) {
            // Skip unavailable grids
            if (!repa::has_grid_type(gt))
                continue;

            // Test consistency of parsing and to_string.
            BOOST_TEST(
                (repa::parse_grid_type(repa::grid_type_to_string(gt)) == gt));

            if (comm.rank() == 0) {
                std::cout << "Checking grid '" << repa::grid_type_to_string(gt)
                          << "'" << std::endl;
            }

            std::unique_ptr<repa::grids::ParallelLCGrid> up = nullptr;
            BOOST_CHECK_NO_THROW(
                up = repa::make_pargrid(gt, comm, box, mings, ep));
            BOOST_TEST(up.get() != nullptr);

            call_test_func(test_func, *this, up.get(), gt);

            if (!repart)
                continue;
            repartition_randomly(up.get());

            call_test_func(test_func, *this, up.get(), gt);

            if (!repart_twice)
                continue;
            repartition_randomly(up.get());

            call_test_func(test_func, *this, up.get(), gt);
        }
    }

    // Calling - with or without grid type:
    // Introspection of the test function if it takes an additional "GridType"
    // parameter.
    using Fno_gt
        = std::function<void(const TEnv &, repa::grids::ParallelLCGrid *)>;
    using Fwith_gt = std::function<void(
        const TEnv &, repa::grids::ParallelLCGrid *, repa::GridType)>;
    static inline void call_test_func(Fno_gt f,
                                      const TEnv &te,
                                      repa::grids::ParallelLCGrid *grid,
                                      repa::GridType gt)
    {
        f(te, grid);
    }

    static inline void call_test_func(Fwith_gt f,
                                      const TEnv &te,
                                      repa::grids::ParallelLCGrid *grid,
                                      repa::GridType gt)
    {
        f(te, grid, gt);
    }
};
} // namespace

inline TEnv default_test_env(repa::ExtraParams ep = repa::ExtraParams{})
{
    return TEnv{ep};
}
