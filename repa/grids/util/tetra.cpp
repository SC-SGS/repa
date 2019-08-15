
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

#define CGAL_DISABLE_ROUNDING_MATH_CHECK
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <boost/iterator/transform_iterator.hpp>

#include "tetra.hpp"

namespace repa {
namespace util {
namespace tetra {

namespace __detail {

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_3<K> Triangulation;
typedef Triangulation::Point Point;
typedef Triangulation::Locate_type Locate_type;

static Point vec2point(const Vec3d &v)
{
    return Point(v[0], v[1], v[2]);
}

struct _Octagon_Impl {
    _Octagon_Impl(const std::array<Vec3d, 8> &vertices);
    const Triangulation T;
};

_Octagon_Impl::_Octagon_Impl(const std::array<Vec3d, 8> &vertices)
    : T(boost::make_transform_iterator(vertices.begin(), vec2point),
        boost::make_transform_iterator(vertices.end(), vec2point))
{
}

static bool contains(const _Octagon_Impl &oi, const Vec3d &v)
{
    Point p = vec2point(v);
    Locate_type lt;
    int li, lj;
    oi.T.locate(p, lt, li, lj);
    // Also accepts corners, edges of the polygon
    return lt <= Triangulation::CELL;
}

} // namespace __detail

// These are declared here because _Octagon_Impl is an imcomplete type in the
// header.
Octagon::Octagon() = default;
Octagon::~Octagon() = default;
Octagon::Octagon(Octagon &&o) = default;

Octagon::Octagon(const std::array<Vec3d, 8> &vertices)
    : oi(std::make_unique<__detail::_Octagon_Impl>(vertices))
{
}

bool Octagon::contains(const Vec3d &p) const
{
    if (!oi)
        throw std::runtime_error("contains() on empty octagon");
    return __detail::contains(*oi, p);
}

bool Octagon::contains(double x, double y, double z) const
{
    return contains({{x, y, z}});
}

void Octagon::operator=(Octagon o)
{
    swap(*this, o);
}

void swap(Octagon &a, Octagon &b)
{
    std::swap(a.oi, b.oi);
}

} // namespace tetra
} // namespace util
} // namespace repa
