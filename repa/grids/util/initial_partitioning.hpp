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

#include "../../pargrid.hpp"
#include "../globox.hpp"
#include "linearize.hpp"
#include "mpi_cart.hpp"
#include "range.hpp"
#include "vec_arith.hpp"
#include <cmath>
#include <stdexcept>

namespace repa {
namespace util {

using namespace repa::grids; // cell index types

/** Enum of possible types of initial partitionings.
 */
enum InitialPartitionType : int {
    LINEAR,
    CARTESIAN1D,
    CARTESIAN2D,
    CARTESIAN3D
};

/** Returns a Cartesian communicator with dims set according to the chosen
 * initial partitioning.
 * If "init_part" does not describe a Cartesian partitioning, returns a
 * null communicator.
 */
boost::mpi::communicator
make_init_part_communicator(const boost::mpi::communicator &comm,
                            InitialPartitionType init_part);

namespace impl {

/** Interface for static partitioning implementations.
 */
template <typename GBox>
struct StaticRankAssigner_impl {
    virtual rank_type rank_of_cell(global_cell_index_type) const = 0;
    virtual ~StaticRankAssigner_impl()
    {
    }
};

/** Assigns a continuous range of cells to each process.
 * The ordering of cells is given by their respective cell index, i.e.
 * by the linearization used.
 */
template <typename GBox>
struct Linear_StaticRankAssigner_impl : public StaticRankAssigner_impl<GBox> {
    Linear_StaticRankAssigner_impl(const GBox &gbox,
                                   const boost::mpi::communicator &comm)
        : _n_cells_per_proc(std::round(
            static_cast<double>(gbox.global_cells().size()) / comm.size())),
          _comm(comm)
    {
    }

    rank_type rank_of_cell(global_cell_index_type g) const override
    {
        return std::min<rank_type>(
            static_cast<rank_type>(g / _n_cells_per_proc), _comm.size() - 1);
    }

private:
    const double _n_cells_per_proc;
    const boost::mpi::communicator &_comm;
};

/** Cartesian partitioning over a box-shaped domain.
 * The shape of the initial partitioning is defined by the dims of "comm".
 * If it is a 1d Cartesian topology, the initial partitioning is 1d, etc.
 */
template <typename GBox>
struct Cart_StaticRankAssigner_impl : public StaticRankAssigner_impl<GBox> {
    Cart_StaticRankAssigner_impl(const GBox &gbox,
                                 const boost::mpi::communicator &comm)
        : _grid_size(gbox.grid_size()),
          _dims(util::mpi_cart_get_dims(comm)),
          _cells_per_proc(_vec3d_div(_grid_size, _dims)),
          _comm(comm)
    {
        assert(comm.has_cartesian_topology());
    }

    rank_type rank_of_cell(global_cell_index_type g) const override
    {
        using namespace repa::util::vector_arithmetic;

        const auto cellidx_3d = util::unlinearize(g, _grid_size);
        const Vec3i procidx3d
            = vec_clamp(_round_vec(cellidx_3d / _cells_per_proc - .5),
                        constant_vec3(0), _dims - 1);
        auto r = util::mpi_cart_rank(_comm, procidx3d);
        assert(r >= 0 && r < _comm.size());
        return r;
    }

private:
    /** Divides two integer vectors as double vectors.
     * This function exists, to make the "using namespace" clause as
     * inner-most as possible.
     */
    static Vec3d _vec3d_div(const Vec3i &a, const Vec3i &b)
    {
        using namespace repa::util::vector_arithmetic;
        return static_cast_vec<Vec3d>(a) / static_cast_vec<Vec3d>(b);
    }

    /** Calls std::round on each component of a double vector and
     * stores the result in an integer vector.
     */
    static Vec3i _round_vec(const Vec3d &a)
    {
        Vec3i res;
        for (int d = 0; d < 3; ++d)
            res[d] = static_cast<int>(std::round(a[d]));
        return res;
    }

    const Vec3i _grid_size;
    const Vec3i _dims;
    const Vec3d _cells_per_proc;
    boost::mpi::communicator
        _comm; // Copy to allow passing an rvalue to the constructor
};

} // namespace impl

/** RAII wrapper class for static rank assigners.
 */
template <typename GBox>
struct StaticRankAssigner {
    StaticRankAssigner(InitialPartitionType pt,
                       const GBox &gbox,
                       const boost::mpi::communicator &comm)
        : _gbox(gbox)
    {
        /** Constructs a StaticRankAssigner_impl from "pt", "gbox" and "comm".
         */
        auto impl_factory = [&]() -> impl::StaticRankAssigner_impl<GBox> * {
            switch (pt) {
            case InitialPartitionType::LINEAR:
                return new impl::Linear_StaticRankAssigner_impl<GBox>(gbox,
                                                                      comm);
            case InitialPartitionType::CARTESIAN1D:
            case InitialPartitionType::CARTESIAN2D:
            case InitialPartitionType::CARTESIAN3D:
                return new impl::Cart_StaticRankAssigner_impl<GBox>(
                    gbox, make_init_part_communicator(comm, pt));
            }
            ensure_not_reached();
        };

        _impl = std::unique_ptr<impl::StaticRankAssigner_impl<GBox>>(
            impl_factory());
    }

    /** Returns the rank of a global cell
     */
    rank_type operator()(global_cell_index_type g) const
    {
        return _impl->rank_of_cell(g);
    }

    /** Returns a range of the initial partitioning (ranks).
     * Ordered the same way as the global cells.
     */
    auto partitioning() const
    {
        return _gbox.global_cells()
               | boost::adaptors::transformed(std::cref(*this));
    }

private:
    std::unique_ptr<impl::StaticRankAssigner_impl<GBox>> _impl;
    const GBox &_gbox;
};

struct UnknownInitialPartitionTypeError : public std::runtime_error {
    UnknownInitialPartitionTypeError(const std::string &s)
        : std::runtime_error("Unknown initial partitioning: " + s)
    {
    }
};

/** Returns the InitialPartitionType with a descriptive string.
 */
util::InitialPartitionType parse_part_type(const std::string &desc);

} // namespace util
} // namespace repa
