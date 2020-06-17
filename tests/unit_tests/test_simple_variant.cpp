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
#define BOOST_TEST_MODULE simple_variant
#include <boost/test/unit_test.hpp>

#include <repa/grids/util/simple_variant.hpp>

using variant_type = repa::util::simple_variant<int, double>;

BOOST_AUTO_TEST_CASE(type_uninitialized)
{
    const variant_type u;
    BOOST_TEST((!u.is_first()));
    BOOST_TEST((!u.is_second()));
    BOOST_TEST((!u.is<int>()));
    BOOST_TEST((!u.is<double>()));
}

BOOST_AUTO_TEST_CASE(type_initialized)
{
    const variant_type u = 1;
    BOOST_TEST(u.is_first());
    BOOST_TEST(u.is<int>());
    BOOST_TEST((!u.is_second()));
    BOOST_TEST((!u.is<double>()));

    const variant_type v = 12.5;
    BOOST_TEST((!v.is_first()));
    BOOST_TEST((!v.is<int>()));
    BOOST_TEST(v.is_second());
    BOOST_TEST(v.is<double>());
}

BOOST_AUTO_TEST_CASE(value)
{
    const variant_type u = 1;

    BOOST_TEST(u.as<int>() == 1);
    BOOST_TEST(static_cast<int>(u) == 1);

    const variant_type v = 12.5;
    BOOST_TEST(v.as<double>() == 12.5);
    BOOST_TEST(static_cast<double>(v) == 12.5);
}

BOOST_AUTO_TEST_CASE(visit)
{
    const variant_type u = 1;

    u.visit([](int i) { BOOST_TEST(i == 1); },
            [](double d) { BOOST_TEST(false); });

    const variant_type v = 12.5;
    v.visit([](int i) { BOOST_TEST(false); },
            [](double d) { BOOST_TEST(d == 12.5); });
}

BOOST_AUTO_TEST_CASE(fmap)
{
    variant_type u = 1;

    double i = u.fmap([](auto x) { return 2 * static_cast<double>(x); });
    BOOST_TEST(i == 2.0);

    variant_type v = 0.5;
    double d = v.fmap([](auto x) -> double { return 2 * x; });
    BOOST_TEST(d == 1.0);
}

BOOST_AUTO_TEST_CASE(as_std_variant)
{
    const variant_type u = 1;

    const auto uv = u.as_std_variant();
    BOOST_TEST(std::holds_alternative<int>(uv));
    BOOST_TEST(std::get<int>(uv) == 1);

    const variant_type v = 12.5;
    const auto vv = v.as_std_variant();
    BOOST_TEST(std::holds_alternative<double>(vv));
    BOOST_TEST(std::get<double>(vv) == 12.5);
}

BOOST_AUTO_TEST_CASE(copy_constructor)
{
    const variant_type u = 1;
    const variant_type v{u};

    BOOST_TEST(u.is_first() == v.is_first());
    BOOST_TEST(u.is_second() == v.is_second());
    BOOST_TEST(u.as<int>() == v.as<int>());
}

BOOST_AUTO_TEST_CASE(assignment)
{
    variant_type u = 1;
    const variant_type v = 12.5;

    u = v;
    BOOST_TEST(u.is_first() == v.is_first());
    BOOST_TEST(u.is_second() == v.is_second());
    BOOST_TEST(u.as<double>() == v.as<double>());
}

BOOST_AUTO_TEST_CASE(comparisons)
{
    const variant_type u = 1;

    const std::vector<int> ints{-1, 0, 1, 2, 3};
    for (const auto i : ints) {
        BOOST_CHECK_EQUAL(u == i, u.as<int>() == i);
        BOOST_CHECK_EQUAL(u != i, u.as<int>() != i);
        BOOST_CHECK_EQUAL(u < i, u.as<int>() < i);
        BOOST_CHECK_EQUAL(u <= i, u.as<int>() <= i);
        BOOST_CHECK_EQUAL(u > i, u.as<int>() > i);
        BOOST_CHECK_EQUAL(u >= i, u.as<int>() >= i);
    }

    const variant_type v = 12.5;

    const std::vector<double> doubles{-10.0, 0.0, 12.4, 12.51, 32.1};
    for (const auto d : doubles) {
        BOOST_CHECK_EQUAL(v == d, v.as<double>() == d);
        BOOST_CHECK_EQUAL(v != d, v.as<double>() != d);
        BOOST_CHECK_EQUAL(v < d, v.as<double>() < d);
        BOOST_CHECK_EQUAL(v <= d, v.as<double>() <= d);
        BOOST_CHECK_EQUAL(v > d, v.as<double>() > d);
        BOOST_CHECK_EQUAL(v >= d, v.as<double>() >= d);
    }

    BOOST_TEST((u != v));
    BOOST_TEST((!(u == v)));
}
