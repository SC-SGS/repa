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

#include <doctest/doctest.h>

#include "../simple_variant.hpp"

using variant_type = repa::util::simple_variant<int, double>;

TEST_CASE("simple_variant type/uninitialized")
{
    const variant_type u;
    CHECK_FALSE(u.is_first());
    CHECK_FALSE(u.is_second());
    CHECK_FALSE(u.is<int>());
    CHECK_FALSE(u.is<double>());
}

TEST_CASE("simple_variant type/initialized")
{
    const variant_type u = 1;
    CHECK(u.is_first());
    CHECK(u.is<int>());
    CHECK_FALSE(u.is_second());
    CHECK_FALSE(u.is<double>());

    const variant_type v = 12.5;
    CHECK_FALSE(v.is_first());
    CHECK_FALSE(v.is<int>());
    CHECK(v.is_second());
    CHECK(v.is<double>());
}

TEST_CASE("simple_variant value")
{
    const variant_type u = 1;

    CHECK(u.as<int>() == 1);
    CHECK(static_cast<int>(u) == 1);

    const variant_type v = 12.5;
    CHECK(v.as<double>() == 12.5);
    CHECK(static_cast<double>(v) == 12.5);
}

TEST_CASE("simple_variant visit")
{
    const variant_type u = 1;

    u.visit([](int i) { CHECK(i == 1); }, [](double d) { CHECK(false); });

    const variant_type v = 12.5;
    v.visit([](int i) { CHECK(false); }, [](double d) { CHECK(d == 12.5); });
}

TEST_CASE("simple_variant fmap")
{
    variant_type u = 1;

    double i = u.fmap([](auto x) { return 2 * static_cast<double>(x); });
    CHECK(i == 2.0);

    variant_type v = 0.5;
    double d = v.fmap([](auto x) -> double { return 2 * x; });
    CHECK(d == 1.0);
}

TEST_CASE("simple_variant std::variant")
{
    const variant_type u = 1;

    const auto uv = u.as_std_variant();
    CHECK(std::holds_alternative<int>(uv));
    CHECK(std::get<int>(uv) == 1);

    const variant_type v = 12.5;
    const auto vv = v.as_std_variant();
    CHECK(std::holds_alternative<double>(vv));
    CHECK(std::get<double>(vv) == 12.5);
}

TEST_CASE("simple_variant copy constructor")
{
    const variant_type u = 1;
    const variant_type v{u};

    CHECK_EQ(u.is_first(), v.is_first());
    CHECK_EQ(u.is_second(), v.is_second());
    CHECK_EQ(u.as<int>(), v.as<int>());
}

TEST_CASE("simple_variant assignment")
{
    variant_type u = 1;
    const variant_type v = 12.5;

    u = v;
    CHECK_EQ(u.is_first(), v.is_first());
    CHECK_EQ(u.is_second(), v.is_second());
    CHECK_EQ(u.as<double>(), v.as<double>());
}

TEST_CASE("simple_variant comparisons")
{
    const variant_type u = 1;

    const std::vector<int> ints{-1, 0, 1, 2, 3};
    for (const auto i : ints) {
        CHECK_EQ(u == i, u.as<int>() == i);
        CHECK_EQ(u != i, u.as<int>() != i);
        CHECK_EQ(u < i, u.as<int>() < i);
        CHECK_EQ(u <= i, u.as<int>() <= i);
        CHECK_EQ(u > i, u.as<int>() > i);
        CHECK_EQ(u >= i, u.as<int>() >= i);
    }

    const variant_type v = 12.5;

    const std::vector<double> doubles{-10.0, 0.0, 12.4, 12.51, 32.1};
    for (const auto d : doubles) {
        CHECK_EQ(v == d, v.as<double>() == d);
        CHECK_EQ(v != d, v.as<double>() != d);
        CHECK_EQ(v < d, v.as<double>() < d);
        CHECK_EQ(v <= d, v.as<double>() <= d);
        CHECK_EQ(v > d, v.as<double>() > d);
        CHECK_EQ(v >= d, v.as<double>() >= d);
    }

    CHECK_FALSE(u == v);
}
