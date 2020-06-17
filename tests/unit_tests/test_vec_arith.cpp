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
#define BOOST_TEST_MODULE vec_arith
#include <boost/test/unit_test.hpp>

#include <repa/common_types.hpp>
#include <repa/grids/util/vec_arith.hpp>

using namespace repa::util::vector_arithmetic;
using repa::Vec3d;
using repa::Vec3i;

BOOST_AUTO_TEST_CASE(binary_ops)
{
    const Vec3i v{7, 9, -21}, w{-1, 2, 17};
    const Vec3i up = v + w, um = v - w, ut = v * w, ud = v / w;
    BOOST_TEST(up == (Vec3i{6, 11, -4}));
    BOOST_TEST(um == (Vec3i{8, 7, -38}));
    BOOST_TEST(ut == (Vec3i{-7, 18, -357}));
    BOOST_TEST(ud == (Vec3i{-7, 4, -1}));
}

BOOST_AUTO_TEST_CASE(binary_ops_literal)
{
    const Vec3i v{7, 9, -21};
    const Vec3i up = v + 2, um = v - 2, ut = v * 2, ud = v / 2;
    BOOST_TEST(up == (Vec3i{9, 11, -19}));
    BOOST_TEST(um == (Vec3i{5, 7, -23}));
    BOOST_TEST(ut == (Vec3i{14, 18, -42}));
    BOOST_TEST(ud == (Vec3i{3, 4, -10}));
}

BOOST_AUTO_TEST_CASE(binary_assignment_ops)
{
    Vec3i v{7, 9, -21}, w{-1, 2, 17};
    v += w;
    BOOST_TEST(v == (Vec3i{6, 11, -4}));
    v -= w;
    BOOST_TEST(v == (Vec3i{7, 9, -21}));
    v *= w;
    BOOST_TEST(v == (Vec3i{-7, 18, -357}));
    v /= w;
    BOOST_TEST(v == (Vec3i{7, 9, -21}));
}

BOOST_AUTO_TEST_CASE(binary_assignment_ops_literal)
{
    Vec3i v{7, 9, -21};
    v += 2;
    BOOST_TEST(v == (Vec3i{9, 11, -19}));
    v -= 2;
    BOOST_TEST(v == (Vec3i{7, 9, -21}));
    v *= 2;
    BOOST_TEST(v == (Vec3i{14, 18, -42}));
    v /= 2;
    BOOST_TEST(v == (Vec3i{7, 9, -21}));
}

BOOST_AUTO_TEST_CASE(shift_operators)
{
    Vec3i v{7, 9, 12}, x{1, 2, 3};
    Vec3i w = v >> 1, y = 1 << x;
    BOOST_TEST(w == (Vec3i{3, 4, 6}));
    BOOST_TEST(y == (Vec3i{2, 4, 8}));
}

BOOST_AUTO_TEST_CASE(unary_operator)
{
    const Vec3i v{7, 9, -21};
    const Vec3i w{-v};
    BOOST_TEST(w == (Vec3i{-7, -9, 21}));
}

BOOST_AUTO_TEST_CASE(expression_capture)
{
    const repa::Vec3i v{1, 1, 1};
    const auto v11 = v + 2;
    const repa::Vec3i v1 = v + 3;
    BOOST_TEST(v1 == (Vec3i{4, 4, 4}));

    const repa::Vec3i v111{v11};
    BOOST_TEST(v111 == (Vec3i{3, 3, 3}));
}

BOOST_AUTO_TEST_CASE(cast_vec)
{
    const Vec3d v{1., 2., 3.};
    const Vec3i w{static_cast_vec<Vec3i>(v)};
    BOOST_TEST(w == (Vec3i{1, 2, 3}));
}

BOOST_AUTO_TEST_CASE(clamp)
{
    const Vec3i v{-2, 1, 10};
    const Vec3i w{vec_clamp(v, constant_vec3(0), constant_vec3(5))};
    BOOST_TEST(w == (Vec3i{0, 1, 5}));
}

BOOST_AUTO_TEST_CASE(wrap)
{
    const Vec3i v{-2, 1, 10};
    const Vec3i w{v % constant_vec3(5)};
    BOOST_TEST(w == (Vec3i{3, 1, 0}));
}

BOOST_AUTO_TEST_CASE(cmp)
{
    const Vec3i v{-2, 1, 10};
    const repa::Vec3<bool> weq{v == -2}, wg{v > -2}, wgeq{v >= -2}, wl{v < -2},
        wleq{v <= -2};
    BOOST_TEST(weq == (repa::Vec3<bool>{true, false, false}));
    BOOST_TEST(wg == (repa::Vec3<bool>{false, true, true}));
    BOOST_TEST(wgeq == (repa::Vec3<bool>{true, true, true}));
    BOOST_TEST(wl == (repa::Vec3<bool>{false, false, false}));
    BOOST_TEST(wleq == (repa::Vec3<bool>{true, false, false}));

    BOOST_TEST(
        all(weq.as_expr() == repa::Vec3<bool>{true, false, false}.as_expr()));
    BOOST_TEST(
        all(wg.as_expr() == repa::Vec3<bool>{false, true, true}.as_expr()));
    BOOST_TEST(
        all(wgeq.as_expr() == repa::Vec3<bool>{true, true, true}.as_expr()));
    BOOST_TEST(
        all(wl.as_expr() == repa::Vec3<bool>{false, false, false}.as_expr()));
    BOOST_TEST(
        all(wleq.as_expr() == repa::Vec3<bool>{true, false, false}.as_expr()));
}
