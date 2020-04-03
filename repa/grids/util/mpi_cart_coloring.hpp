/**
 * Copyright 2017-2020 Steffen Hirschmann
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

#include "mpi_cart.hpp"

namespace repa {
namespace util {

struct independent_process_sets {
    independent_process_sets(const boost::mpi::communicator &comm_cart)
        : _coords(util::mpi_cart_get_coords(comm_cart))
    {
        assert(comm_cart.has_cartesian_topology());
    }

    template <typename Func>
    independent_process_sets &for_each(Func &&f)
    {
        _round_func = f;
        return *this;
    }

    template <typename Func>
    independent_process_sets &for_all_after_each_round(Func &&f)
    {
        _after_func = f;
        return *this;
    }

    void operator()()
    {
        for (int color = 0; color < 8; color++) {
            if (_coords[0] % 2 == !!(color & 0x1)
                && _coords[1] % 2 == !!(color & 0x2)
                && _coords[2] % 2 == !!(color & 0x4)) {
                _round_func();
            }
            _after_func();
        }
    }

private:
    const Vec3i _coords;
    std::function<void(void)> _round_func, _after_func;
};

} // namespace util
} // namespace repa