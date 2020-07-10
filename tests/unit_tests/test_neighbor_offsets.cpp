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
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE neighbor_offsets
#include <boost/test/unit_test.hpp>

#include <algorithm>

#include <repa/grids/util/neighbor_offsets.hpp>

BOOST_AUTO_TEST_CASE(uniqueness)
{
    // Create copy
    auto no = repa::util::NeighborOffsets3D::raw;
    std::sort(no.begin(), no.end(), [](const auto &a, const auto &b) {
        // Use lexicographical sorting of std::array.
        return a.as_array() < b.as_array();
    });
    auto it = std::unique(no.begin(), no.end());
    BOOST_TEST(it == no.end());
}

BOOST_AUTO_TEST_CASE(half_shell)
{
    const auto &no = repa::util::NeighborOffsets3D::raw;
    const auto hs_begin = no.begin() + 1;
    const auto hs_end = no.begin() + 14;

    // Additive inverse not allowed to be included in the half-shell
    for (int i = 1; i <= 13; ++i) {
        repa::Vec3i v = no[i];
        for (int d = 0; d < 3; ++d)
            v[d] = -v[d];
        BOOST_TEST(std::find(hs_begin, hs_end, v) == hs_end);
    }
}

BOOST_AUTO_TEST_CASE(order)
{
    BOOST_TEST(
        (repa::util::neighbor_type(0) == repa::util::NeighborCellType::SELF));
    for (int i = 1; i <= 13; ++i)
        BOOST_TEST((repa::util::neighbor_type(i)
                    == repa::util::NeighborCellType::HALF_SHELL));
    for (int i = 14; i < 27; ++i)
        BOOST_TEST((repa::util::neighbor_type(i)
                    == repa::util::NeighborCellType::FULL_SHELL));
}