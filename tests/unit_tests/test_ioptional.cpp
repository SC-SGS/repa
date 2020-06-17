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
#define BOOST_TEST_MODULE ioptional
#include <boost/test/unit_test.hpp>

#include <map>
#include <string>
#include <unordered_map>

#include <repa/grids/util/ioptional.hpp>

using optional_int = repa::util::ioptional<repa::rank_type>;

BOOST_AUTO_TEST_CASE(empty)
{
    auto io = optional_int{};
    BOOST_TEST((!io.has_value()));
    io = 12;
    BOOST_TEST((io.has_value()));

    auto io2 = optional_int{12};
    BOOST_TEST((io2.has_value()));
}

BOOST_AUTO_TEST_CASE(value)
{
    auto io = optional_int{};
    BOOST_TEST((io.value_or(-255) == -255));

    BOOST_CHECK_THROW(io.value_or_throw<std::runtime_error>("test"),
                      std::runtime_error);

    io = 12;
    BOOST_TEST((io.value() == 12));
    BOOST_TEST((io.value_or(-255) == 12));
    BOOST_TEST((io.ensure_value() == 12));
    BOOST_TEST((io.value_or_throw<std::runtime_error>("test") == 12));
}

BOOST_AUTO_TEST_CASE(operator_bool)
{
    auto io = optional_int{};
    BOOST_TEST((!static_cast<bool>(io)));
    io = 12;
    BOOST_TEST(static_cast<bool>(io));
}

BOOST_AUTO_TEST_CASE(deref)
{
    auto io = optional_int{12};
    BOOST_TEST((*io == 12));
}

BOOST_AUTO_TEST_CASE(reset)
{
    auto io = optional_int{12};
    BOOST_TEST((io.has_value()));
    io.reset();
    BOOST_TEST((!io.has_value()));
}

BOOST_AUTO_TEST_CASE(cmp)
{
    auto io = optional_int{};
    auto io2 = optional_int{12};
    BOOST_TEST((io != io2));

    io = io2;
    BOOST_TEST((io == io2));
}