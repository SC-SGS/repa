
/**
 * Copyright 2017-2020 Steffen Hirschmann
 * Copyright 2020 Benjamin Vier
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
#include <cmath>

namespace repa {
namespace util {
namespace tetra {

using Vec3i64 = Vec3<int64_t>;

// Anonymous namespace for internal linkage
namespace {

using namespace vector_arithmetic;

Vec3i64 integerize(const Vec3d &v)
{
    return static_cast_vec<Vec3i64>(v * static_cast<double>(precision));
}

std::array<Vec3i64, 8> integerizedArray(const std::array<Vec3d, 8> &vertices)
{
    std::array<Vec3i64, 8> intVert;
    for (int i = 0; i < 8; i++) {
        intVert[i] = integerize(vertices[i]);
    }
    return intVert;
}

struct Plane {
    Vec3i64 normVector;
    int64_t heightOfPlane;

    Plane(){};

    Plane(const std::array<Vec3i64, 3> &vecs)
    {
        normVector = cross(vecs[0] - vecs[2], vecs[1] - vecs[0]);
        heightOfPlane = dot(normVector, vecs[0]);
    }

    bool isAboveOrEqual(Vec3i64 point) const noexcept
    {
        return dot(point, normVector) >= heightOfPlane;
    }

    bool isAbove(Vec3i64 point) const noexcept
    {
        return dot(point, normVector) > heightOfPlane;
    }

    bool is2RcAbove(Vec3i64 point, double twoRc) const noexcept
    {
        int height2Rc = static_cast<int>(
            ceil(twoRc * static_cast<double>(sumOfAbs(normVector))));
        return dot(point, normVector) > heightOfPlane + height2Rc;
    }
};

} // namespace

struct _Octagon_Impl {
private:
    static const std::array<int, 6> cornerOrder;
    Plane tetras[6][4];
    double min_height;
    bool isValid;

public:
    _Octagon_Impl() = delete;

    _Octagon_Impl(const std::array<Vec3i64, 8> &corners, double max_cutoff)
        : min_height(2. * std::sqrt(3.) * max_cutoff), isValid(true)
    {
        Vec3i64 start = corners[0];
        Vec3i64 end = corners[7];
        Vec3i64 last = corners[5];
        for (int i = 0; i < 6; i++) {
            Vec3i64 next = corners[cornerOrder[i]];
            addTetra(i, {start, end, next, last});
            last = next;
        }
    }

    void addTetra(int index, const std::array<Vec3i64, 4> &corners)
    {
        tetras[index][0] = Plane({corners[0], corners[1], corners[2]});
        tetras[index][1] = Plane({corners[0], corners[2], corners[3]});
        tetras[index][2] = Plane({corners[0], corners[3], corners[1]});
        tetras[index][3] = Plane({corners[1], corners[3], corners[2]});
        if (has_validity_check() && isValid) {
            isValid = isValid
                      && tetras[index][0].is2RcAbove(corners[3], min_height);
            isValid = isValid
                      && tetras[index][1].is2RcAbove(corners[1], min_height);
            isValid = isValid
                      && tetras[index][2].is2RcAbove(corners[2], min_height);
            isValid = isValid
                      && tetras[index][3].is2RcAbove(corners[0], min_height);
        }
    }

    bool contains(Vec3i64 point) const noexcept
    {
        // Iterate over all tetrahedrons of the domain
        for (int tetra = 0; tetra < 6; tetra++) {
            // Points which are exactly on this plane of the tetrahedron
            // are not accepted to avoid that they are assigned to two domains.
            bool tetraContainsP = tetras[tetra][3].isAbove(point);
            // Iterate over all other sides
            for (int plane = 0; plane < 3; plane++) {
                if (!tetraContainsP) {
                    break;
                }
                tetraContainsP = tetras[tetra][plane].isAboveOrEqual(point);
            }
            // The point is accepted when it lies within a tetrahedron of the
            // domain
            if (tetraContainsP) {
                return true;
            }
        }
        return false;
    }

    bool has_validity_check() const noexcept
    {
        return min_height > 0.0;
    }

    bool is_valid() const
    {
        return isValid;
    }
};

const std::array<int, 6> _Octagon_Impl::cornerOrder = {{1, 3, 2, 6, 4, 5}};

// These are declared here because _Octagon_Impl is an imcomplete type in the
// header.
Octagon::Octagon() = default;
Octagon::~Octagon() = default;
Octagon::Octagon(Octagon &&o) = default;

Octagon::Octagon(const std::array<Vec3d, 8> &vertices)
    : oi(std::make_unique<_Octagon_Impl>(integerizedArray(vertices), 0.0))
{
}

Octagon::Octagon(const std::array<Vec3d, 8> &vertices, double max_cs)
    : oi(std::make_unique<_Octagon_Impl>(integerizedArray(vertices), max_cs))
{
}

bool Octagon::is_valid() const
{
    if (!oi)
        throw std::runtime_error("is_valid() on empty octagon");
    if (!oi->has_validity_check())
        throw std::runtime_error(
            "The validity of this octagon wasn't checked on initialization");
    return oi->is_valid();
}

bool Octagon::contains(const Vec3d &p) const
{
    if (!oi)
        throw std::runtime_error("contains() on empty octagon");
    return oi->contains(integerize(p));
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
