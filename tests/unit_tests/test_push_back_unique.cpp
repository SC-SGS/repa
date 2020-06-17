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
#define BOOST_TEST_MODULE push_back_unique
#include <boost/test/unit_test.hpp>

#include "repa/grids/util/push_back_unique.hpp"

BOOST_AUTO_TEST_CASE(POD_data)
{
    std::vector<int> v{1, 2, 3, 4};
    const auto v_original = v;

    {
        repa::util::push_back_unique(v, 1);
        repa::util::push_back_unique(v, 4);
        BOOST_TEST(v == v_original);
    }

    {
        repa::util::push_back_unique(v, 5);
        repa::util::push_back_unique(v, 6);
        BOOST_TEST(v.size() == v_original.size() + 2);

        std::vector<int> comp = v_original;
        comp.push_back(5);
        comp.push_back(6);
        BOOST_TEST(v == comp);
    }
}

struct Empty {
};

bool operator==(const Empty &a, const Empty &b)
{
    (void)a;
    (void)b;
    return true;
}

BOOST_AUTO_TEST_CASE(empty_custom_struct)
{
    std::vector<Empty> v;

    repa::util::push_back_unique(v, Empty{});
    BOOST_TEST(v.size() == 1);

    repa::util::push_back_unique(v, Empty{});
    BOOST_TEST(v.size() == 1);
}
