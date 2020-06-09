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

#include "common_types.hpp"
#include "../vec_arith.hpp"

using namespace repa::util::vector_arithmetic;
using repa::Vec3d;
using repa::Vec3i;

TEST_CASE("vec_arith binary ops")
{
    const Vec3i v{7, 9, -21}, w{-1, 2, 17};
    const Vec3i up = v + w, um = v - w, ut = v * w, ud = v / w;
    CHECK(up == Vec3i{6, 11, -4});
    CHECK(um == Vec3i{8, 7, -38});
    CHECK(ut == Vec3i{-7, 18, -357});
    CHECK(ud == Vec3i{-7, 4, -1});
}

TEST_CASE("vec_arith binary ops literal")
{
    const Vec3i v{7, 9, -21};
    const Vec3i up = v + 2, um = v - 2, ut = v * 2, ud = v / 2;
    CHECK(up == Vec3i{9, 11, -19});
    CHECK(um == Vec3i{5, 7, -23});
    CHECK(ut == Vec3i{14, 18, -42});
    CHECK(ud == Vec3i{3, 4, -10});
}

TEST_CASE("vec_arith binary assignment ops")
{
    Vec3i v{7, 9, -21}, w{-1, 2, 17};
    v += w;
    CHECK(v == Vec3i{6, 11, -4});
    v -= w;
    CHECK(v == Vec3i{7, 9, -21});
    v *= w;
    CHECK(v == Vec3i{-7, 18, -357});
    v /= w;
    CHECK(v == Vec3i{7, 9, -21});
}

TEST_CASE("vec_arith binary assignment ops literal")
{
    Vec3i v{7, 9, -21};
    v += 2;
    CHECK(v == Vec3i{9, 11, -19});
    v -= 2;
    CHECK(v == Vec3i{7, 9, -21});
    v *= 2;
    CHECK(v == Vec3i{14, 18, -42});
    v /= 2;
    CHECK(v == Vec3i{7, 9, -21});
}

TEST_CASE("vec_arith shift operators")
{
    Vec3i v{7, 9, 12}, x{1, 2, 3};
    Vec3i w = v >> 1, y = 1 << x;
    CHECK(w == Vec3i{3, 4, 6});
    CHECK(y == Vec3i{2, 4, 8});
}

TEST_CASE("vec_arith unary operator")
{
    const Vec3i v{7, 9, -21};
    const Vec3i w{-v};
    CHECK(w == Vec3i{-7, -9, 21});
}

TEST_CASE("vec_arith expression capture")
{
    const repa::Vec3i v{1, 1, 1};
    const auto v11 = v + 2;
    const repa::Vec3i v1 = v + 3;
    CHECK(v1 == Vec3i{4, 4, 4});

    const repa::Vec3i v111{v11};
    CHECK(v111 == Vec3i{3, 3, 3});
}

TEST_CASE("vec_arith static_cast_vec")
{
    const Vec3d v{1., 2., 3.};
    const Vec3i w{static_cast_vec<Vec3i>(v)};
    CHECK(w == Vec3i{1, 2, 3});
}

TEST_CASE("vec_arith clamp")
{
    const Vec3i v{-2, 1, 10};
    const Vec3i w{vec_clamp(v, constant_vec3(0), constant_vec3(5))};
    CHECK(w == Vec3i{0, 1, 5});
}

TEST_CASE("vec_arith wrap")
{
    const Vec3i v{-2, 1, 10};
    const Vec3i w{v % constant_vec3(5)};
    CHECK(w == Vec3i{3, 1, 0});
}

TEST_CASE("vec_arith cmp")
{
    const Vec3i v{-2, 1, 10};
    const repa::Vec3<bool> weq{v == -2}, wg{v > -2}, wgeq{v >= -2}, wl{v < -2},
        wleq{v <= -2};
    CHECK(weq == repa::Vec3<bool>{true, false, false});
    CHECK(wg == repa::Vec3<bool>{false, true, true});
    CHECK(wgeq == repa::Vec3<bool>{true, true, true});
    CHECK(wl == repa::Vec3<bool>{false, false, false});
    CHECK(wleq == repa::Vec3<bool>{true, false, false});

    CHECK(all(weq.as_expr() == repa::Vec3<bool>{true, false, false}.as_expr()));
    CHECK(all(wg.as_expr() == repa::Vec3<bool>{false, true, true}.as_expr()));
    CHECK(all(wgeq.as_expr() == repa::Vec3<bool>{true, true, true}.as_expr()));
    CHECK(all(wl.as_expr() == repa::Vec3<bool>{false, false, false}.as_expr()));
    CHECK(
        all(wleq.as_expr() == repa::Vec3<bool>{true, false, false}.as_expr()));
}
