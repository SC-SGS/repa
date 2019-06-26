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
#include <boost/mpi/environment.hpp>
#include <random>

struct TestEnv {
    boost::mpi::environment env;
    boost::mpi::communicator comm;
    repa::Vec3d box;
    double mings;

    TestEnv(const repa::Vec3d &box, double min_grid_size)
        : box(box), mings(min_grid_size)
    {
    }

    const boost::mpi::communicator &get_comm()
    {
        return comm;
    }

protected:
    template <typename Func, typename PPFunc, typename... Args>
    void __run_for_all_grid_type_impl(repa::GridType gt,
                                      Func testf,
                                      PPFunc preprocess,
                                      Args... args)
    {
        if (comm.rank() == 0) {
            std::cout << "Checking grid '" << repa::grid_type_to_string(gt)
                      << "'" << std::endl;
        }

        std::unique_ptr<repa::grids::ParallelLCGrid> up = nullptr;
        BOOST_CHECK_NO_THROW(up = repa::make_pargrid(gt, comm, box, mings));
        BOOST_TEST(up.get() != nullptr);

        preprocess(up.get());
        testf(up.get(), args...);
    }

    template <typename Func, typename PPFunc, typename... Args>
    void __run_for_all_all_grid_types_impl(Func testf,
                                           PPFunc preprocess,
                                           Args... args)
    {
        for (const auto gt : repa::supported_grid_types()) {
            __run_for_all_grid_type_impl(gt, testf, preprocess, args...);
        }
    }
};

struct StaticTestEnv : public TestEnv {
    StaticTestEnv(const repa::Vec3d &box, double min_grid_size)
        : TestEnv(box, min_grid_size)
    {
    }

    template <typename Func, typename... Args>
    void run_for_all_grid_types(Func testf, Args... args)
    {
        auto no_preprocessing = [](repa::grids::ParallelLCGrid *) {};
        __run_for_all_all_grid_types_impl(testf, no_preprocessing, args...);
    }
};

static void repartition_randomly(repa::grids::ParallelLCGrid *grid)
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
    auto random_dbl2 = [&random_dbl](int, int) { return random_dbl(); };

    BOOST_CHECK_NO_THROW(
        grid->repartition(random_cellweights, random_dbl2, noop));
}

struct RepartTestEnv : public TestEnv {
    RepartTestEnv(const repa::Vec3d &box, double min_grid_size)
        : TestEnv(box, min_grid_size)
    {
    }

    template <typename Func, typename... Args>
    void run_for_all_grid_types(Func testf, Args... args)
    {
        auto no_preprocessing = [](repa::grids::ParallelLCGrid *) {};
        __run_for_all_all_grid_types_impl(testf, no_preprocessing, args...);
        __run_for_all_all_grid_types_impl(testf, repartition_randomly, args...);
    }
};