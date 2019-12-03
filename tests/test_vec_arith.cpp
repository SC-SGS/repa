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

/**
 * Checks vec_arith.hpp
 */

#define BOOST_TEST_MODULE vec_arith

#include <boost/test/included/unit_test.hpp>
#include <repa/common_types.hpp>
#include <repa/grids/util/vec_arith.hpp>

using namespace repa::util::vector_arithmetic;
using repa::Vec3d;
using repa::Vec3i;

void check_binary_ops()
{
    const Vec3i v{7, 9, -21}, w{-1, 2, 17};
    const Vec3i up = v + w, um = v - w, ut = v * w, ud = v / w;
    BOOST_CHECK((up == Vec3i{6, 11, -4}));
    BOOST_CHECK((um == Vec3i{8, 7, -38}));
    BOOST_CHECK((ut == Vec3i{-7, 18, -357}));
    BOOST_CHECK((ud == Vec3i{-7, 4, -1}));
}

void check_binary_ops_literal()
{
    const Vec3i v{7, 9, -21};
    const Vec3i up = v + 2, um = v - 2, ut = v * 2, ud = v / 2;
    BOOST_CHECK((up == Vec3i{9, 11, -19}));
    BOOST_CHECK((um == Vec3i{5, 7, -23}));
    BOOST_CHECK((ut == Vec3i{14, 18, -42}));
    BOOST_CHECK((ud == Vec3i{3, 4, -10}));
}

void check_binary_assignment_ops()
{
    Vec3i v{7, 9, -21}, w{-1, 2, 17};
    v += w;
    BOOST_CHECK((v == Vec3i{6, 11, -4}));
    v -= w;
    BOOST_CHECK((v == Vec3i{7, 9, -21}));
    v *= w;
    BOOST_CHECK((v == Vec3i{-7, 18, -357}));
    v /= w;
    BOOST_CHECK((v == Vec3i{7, 9, -21}));
}

void check_binary_assignment_ops_literal()
{
    Vec3i v{7, 9, -21};
    v += 2;
    BOOST_CHECK((v == Vec3i{9, 11, -19}));
    v -= 2;
    BOOST_CHECK((v == Vec3i{7, 9, -21}));
    v *= 2;
    BOOST_CHECK((v == Vec3i{14, 18, -42}));
    v /= 2;
    BOOST_CHECK((v == Vec3i{7, 9, -21}));
}

void check_unary()
{
    const Vec3i v{7, 9, -21};
    const Vec3i w{-v};
    BOOST_CHECK((w == Vec3i{-7, -9, 21}));
}

void check_expression_capture()
{
    const repa::Vec3i v{1, 1, 1};
    const auto v11 = v + 2;
    const repa::Vec3i v1 = v + 3;
    BOOST_CHECK((v1 == Vec3i{4, 4, 4}));

    const repa::Vec3i v111{v11};
    BOOST_CHECK((v111 == Vec3i{3, 3, 3}));
}

void check_cast()
{
    const Vec3d v{1., 2., 3.};
    const Vec3i w{static_cast_vec<Vec3i>(v)};
    BOOST_CHECK((w == Vec3i{1, 2, 3}));
}

void check_clamp()
{
    const Vec3i v{-2, 1, 10};
    const Vec3i w{vec_clamp(v, constant_vec3(0), constant_vec3(5))};
    BOOST_CHECK((w == Vec3i{0, 1, 5}));
}

void check_wrap()
{
    const Vec3i v{-2, 1, 10};
    const Vec3i w{v % constant_vec3(5)};
    BOOST_CHECK((w == Vec3i{3, 1, 0}));
}

void check_comparison()
{
    const Vec3i v{-2, 1, 10};
    const repa::Vec3<bool> weq{v == -2}, wg{v > -2}, wgeq{v >= -2}, wl{v < -2},
        wleq{v <= -2};
    BOOST_CHECK((weq == repa::Vec3<bool>{true, false, false}));
    BOOST_CHECK((wg == repa::Vec3<bool>{false, true, true}));
    BOOST_CHECK((wgeq == repa::Vec3<bool>{true, true, true}));
    BOOST_CHECK((wl == repa::Vec3<bool>{false, false, false}));
    BOOST_CHECK((wleq == repa::Vec3<bool>{true, false, false}));

    BOOST_CHECK(
        all(weq.as_expr() == repa::Vec3<bool>{true, false, false}.as_expr()));
    BOOST_CHECK(
        all(wg.as_expr() == repa::Vec3<bool>{false, true, true}.as_expr()));
    BOOST_CHECK(
        all(wgeq.as_expr() == repa::Vec3<bool>{true, true, true}.as_expr()));
    BOOST_CHECK(
        all(wl.as_expr() == repa::Vec3<bool>{false, false, false}.as_expr()));
    BOOST_CHECK(
        all(wleq.as_expr() == repa::Vec3<bool>{true, false, false}.as_expr()));
}

BOOST_AUTO_TEST_CASE(test_vec_arith)
{
    check_binary_ops();
    check_binary_ops_literal();
    check_binary_assignment_ops();
    check_binary_assignment_ops_literal();
    check_unary();
    check_expression_capture();
    check_cast();
    check_clamp();
    check_wrap();
    check_comparison();
}
