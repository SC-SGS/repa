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
#define BOOST_TEST_MODULE integral_range
#include <boost/test/unit_test.hpp>

#include <repa/common_types.hpp>

BOOST_AUTO_TEST_CASE(domain)
{
    for (repa::fs_neighidx::value_type i = 0; i <= 26; ++i) {
        repa::fs_neighidx ni{i};
        BOOST_TEST(static_cast<repa::fs_neighidx::value_type>(ni) == i);
    }

#ifndef NDEBUG
    repa::fs_neighidx ni = 0;
    BOOST_CHECK_THROW(ni = 27, std::domain_error);
    BOOST_CHECK_THROW(ni = -1, std::domain_error);
#endif
}
