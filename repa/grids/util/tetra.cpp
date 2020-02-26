
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

#include "tetra.hpp"
#include "vec_arith.hpp"

namespace repa {
namespace util {
namespace tetra {

namespace __detail {

Vec3i integerize(Vec3d v)
{
    using namespace vector_arithmetic;
    v *= precision;
    return static_cast_vec<Vec3i>(v);
}

std::array<Vec3i, 8> integerizedArray(std::array<Vec3d, 8> vertices)
{
    std::array<Vec3i, 8> intVert = {};
    for (int i = 0; i < 8; i++) {
        intVert[i] = integerize(vertices[i]);
    }
    return intVert;
}

struct Plane {
    Vec3i normVector;
    int64_t heightOfPlane;
    Plane(){};
    Plane(std::array<Vec3i, 3> vecs);
    bool isAboveOrEqual(Vec3i v);
    bool isAbove(Vec3i v);
};

Plane::Plane(std::array<Vec3i, 3> vecs)
{
    using namespace vector_arithmetic;
    Vec3i v1 = vecs[1] - vecs[0];
    Vec3i v2 = vecs[0] - vecs[2];
    this->normVector = cross(v2, v1);
    this->heightOfPlane = dot(this->normVector, vecs[0]);
}

bool Plane::isAboveOrEqual(Vec3i point)
{
    using vector_arithmetic::dot;
    return dot(point, this->normVector) >= this->heightOfPlane;
}

bool Plane::isAbove(Vec3i point)
{
    using vector_arithmetic::dot;
    return dot(point, this->normVector) > this->heightOfPlane;
}

struct _Octagon_Impl {
    int cornerOrder[6] = {1, 3, 2, 6, 4, 5};
    Plane tetras[6][4];

public:
    _Octagon_Impl(std::array<Vec3i, 8> corners);
    void addTetra(int index, std::array<Vec3i, 4> tetraCorners);
};
_Octagon_Impl::_Octagon_Impl(std::array<Vec3i, 8> corners)
{
    Vec3i start = corners[0];
    Vec3i end = corners[7];
    Vec3i last = corners[5];
    for (int i = 0; i < 6; i++) {
        Vec3i next = corners[this->cornerOrder[i]];
        this->addTetra(i, {start, end, last, next});
        last = next;
    }
}

void _Octagon_Impl::addTetra(int index, std::array<Vec3i, 4> corners)
{
    this->tetras[index][0] = Plane({corners[0], corners[1], corners[2]});
    this->tetras[index][1] = Plane({corners[0], corners[2], corners[3]});
    this->tetras[index][2] = Plane({corners[0], corners[3], corners[1]});
    this->tetras[index][3] = Plane({corners[1], corners[3], corners[2]});
}

bool contains(_Octagon_Impl oi, Vec3i point)
{
    // Iterate over all tetrahedrons of the domain
    for (int tetra = 0; tetra < 6; tetra++) {
        // Points which are exactly on this plane of the tetrahedron
        // are not accepted to avoid that they are assigned to two domains.
        bool tetraContainsP = oi.tetras[tetra][3].isAbove(point);
        // Iterate over all other sides
        for (int plane = 0; plane < 3; plane++) {
            if (!tetraContainsP) {
                break;
            }
            tetraContainsP = oi.tetras[tetra][plane].isAboveOrEqual(point);
        }
        // The point is accepted when it lies within a tetrahedron of the domain
        if (tetraContainsP) {
            return true;
        }
    }
    return false;
}

} // namespace __detail

// These are declared here because _Octagon_Impl is an imcomplete type in the
// header.
Octagon::Octagon() = default;
Octagon::~Octagon() = default;
Octagon::Octagon(Octagon &&o) = default;

Octagon::Octagon(const std::array<Vec3d, 8> &vertices)
    : oi(std::make_unique<__detail::_Octagon_Impl>(
        __detail::integerizedArray(vertices)))
{
}

bool Octagon::contains(const Vec3d &p) const
{
    if (!oi)
        throw std::runtime_error("contains() on empty octagon");
    return __detail::contains(*oi, __detail::integerize(p));
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
