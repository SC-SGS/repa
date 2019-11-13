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
#include "testenv.hpp"
#include <boost/test/unit_test.hpp>
#include <optional>
#include <random>

testenv::TEnv default_test_env(repa::ExtraParams ep)
{
    return testenv::TEnv{ep};
}

static std::vector<double> get_random_vec(size_t n)
{
    std::random_device rd;
    std::mt19937 gen{rd()};
    std::uniform_real_distribution<> dis{1., 10.};

    std::vector<double> v;
    v.reserve(n);
    std::generate_n(std::back_inserter(v), n,
                    [&dis, &gen]() { return dis(gen); });

    return v;
}

static void repartition_helper(repa::grids::ParallelLCGrid *grid,
                               testenv::MetricFunc *f)
{
    auto noop = []() {};
    auto ccm = [](int, int) { return 1.0; };
    std::vector<double> metric_values
        = (f ? *f : get_random_vec)(static_cast<size_t>(grid->n_local_cells()));
    auto metric = [&metric_values]() { return metric_values; };

    BOOST_CHECK_NO_THROW(grid->repartition(metric, ccm, noop));
}

namespace testenv {

TEnv::TEnv(const boost::mpi::communicator &comm,
           repa::Vec3d box,
           double mings,
           repa::ExtraParams ep)
    : comm(comm, boost::mpi::comm_duplicate),
      box(std::move(box)),
      mings(mings),
      ep(ep)
{
}

TEnv::TEnv(repa::ExtraParams ep) : mings(1.0), ep(ep)
{
    // Devise some appropriately sized grid suitable for all methods.
    repa::Vec3i dims{0, 0, 0};
    MPI_Dims_create(comm.size(), 3, dims.data());
    for (size_t i = 0; i < box.size(); ++i)
        box[i] = mings * 5 * (dims[i] + 1);
}

TEnv &TEnv::with_repart()
{
    repart = true;
    repart_twice = false;
    return *this;
}

TEnv &TEnv::with_repart(MetricFunc &f)
{
    get_metric = &f;
    return with_repart();
}

TEnv &TEnv::with_repart_twice()
{
    repart = true;
    repart_twice = true;
    return *this;
}

TEnv &TEnv::without_repart()
{
    repart = false;
    repart_twice = false;
    return *this;
}

TEnv &TEnv::all_grids()
{
    grids = repa::supported_grid_types();
    return *this;
}

TEnv &TEnv::only(std::set<repa::GridType> s)
{
    grids.clear();
    auto supported = repa::supported_grid_types();
    for (const auto gt : s) {
        if (supported.find(gt) != std::end(supported))
            grids.insert(gt);
    }
    return *this;
}

TEnv &TEnv::exclude(std::set<repa::GridType> s)
{
    for (const auto gt : s)
        grids.erase(gt);
    return *this;
}

void TEnv::run_impl(TestFunc test_func)
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
        BOOST_CHECK_NO_THROW(up = repa::make_pargrid(gt, comm, box, mings, ep));
        BOOST_TEST(up.get() != nullptr);

        test_func(*this, up.get(), gt);

        if (!repart)
            continue;
        repartition_helper(up.get(), get_metric);

        test_func(*this, up.get(), gt);

        if (!repart_twice)
            continue;
        repartition_helper(up.get(), get_metric);

        test_func(*this, up.get(), gt);
    }
}

void TEnv::run(Fno_gt test_func)
{
    run_impl([test_func](auto te, auto grid, auto) { test_func(te, grid); });
}
void TEnv::run(Fwith_gt test_func)
{
    run_impl(
        [test_func](auto te, auto grid, auto gt) { test_func(te, grid, gt); });
}
} // namespace testenv
