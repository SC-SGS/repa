
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
#include <iostream>

namespace repa {
namespace util {
namespace tetra {

int16_t precision = 10;
Vec3i box_size{0, 0, 0};

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

std::pair<Vec3i, Vec3i> min_max_of_dim(const std::array<Vec3i64, 8> &vertices)
{
    Vec3i min = box_size;
    Vec3i max{0, 0, 0};
    for (int d = 0; d < 3; d++) {
        for (Vec3i64 vertex : vertices) {
            int v_d = vertex[d];
            min[d] = std::min(v_d, min[d]);
            max[d] = std::max(v_d, max[d]);
        }
    }
    return std::make_pair(min, max);
}

struct Plane {
    Vec3i64 normVector;
    int64_t heightOfPlane;

    Plane()
    {
    }

    Plane(const std::array<Vec3i64, 3> &vecs)
        : normVector(cross(vecs[0] - vecs[2], vecs[1] - vecs[0])),
          heightOfPlane(dot(normVector, vecs[0]))
    {
    }

    bool isAboveOrEqual(Vec3i64 point) const noexcept
    {
        return dot(point, normVector) >= heightOfPlane;
    }

    bool isAbove(Vec3i64 point) const noexcept
    {
        return dot(point, normVector) > heightOfPlane;
    }

    bool isXAbove(Vec3i64 point, double dist) const noexcept
    {
        int64_t int_dist = static_cast<int64_t>(
            ceil(dist * precision * static_cast<double>(sumOfAbs(normVector))));
        return dot(point, normVector) > heightOfPlane + int_dist;
    }
};

std::array<Plane, 4>
planes_of_tetrahedron(const std::array<Vec3i64, 4> &corners)
{
    return {Plane({corners[0], corners[1], corners[2]}),
            Plane({corners[0], corners[2], corners[3]}),
            Plane({corners[0], corners[3], corners[1]}),
            Plane({corners[1], corners[3], corners[2]})};
}

} // namespace

struct _Octagon_Impl {
private:
    static const std::array<int, 6> cornerOrder;
    std::array<std::array<Plane, 4>, 6> tetrahedron_planes;
    double min_height;
    bool isValid;
    Vec3<bool> periodic;
    Vec3i min_dim;

public:
    _Octagon_Impl() = delete;

    _Octagon_Impl(const std::array<Vec3i64, 8> &corners, double max_cutoff)
        : min_height(2. * std::sqrt(3.) * max_cutoff), isValid(true)
    {
        periodic = {false, false, false};
        if (util::vector_arithmetic::all(box_size > 0)) {
            check_periodity(corners);
        }
        Vec3i64 start = corners[0];
        Vec3i64 end = corners[7];
        Vec3i64 last = corners[5];
        for (int i = 0; i < 6; i++) {
            Vec3i64 next = corners[cornerOrder[i]];
            tetrahedron_planes[i] = generate_tetra({start, end, next, last});
            last = next;
        }
    }

    void check_periodity(const std::array<Vec3i64, 8> corners)
    {
        auto minmax = min_max_of_dim(corners);
        Vec3i min = minmax.first;
        Vec3i max = minmax.second;
        for (int d = 0; d < 3; d++) {
            if (max[d] - min[d] > box_size[d]) {
                std::cerr << "A domain with a lenght greater than the box "
                             "itself is technically possible but not accepted.";
            }
            if (max[d] > box_size[d]) { // box_size -1?
                periodic[d] = true;
            }
            if (min[d] < 0) {
                for (Vec3i64 vec : corners) {
                    vec[d] += box_size[d];
                }
                min[d] += box_size[d];
                periodic[d] = true;
            }
        }
        min_dim = min;
    }

    /** Returns 4 planes that represent the faces of a tetrahedron.
     * Additionally, updates isValid.
     */
    inline std::array<Plane, 4>
    generate_tetra(const std::array<Vec3i64, 4> &corners)
    {
        auto cur_tetra = planes_of_tetrahedron(corners);
        if (has_validity_check())
            isValid = isValid && cur_tetra[0].isXAbove(corners[3], min_height)
                      && cur_tetra[1].isXAbove(corners[1], min_height)
                      && cur_tetra[2].isXAbove(corners[2], min_height)
                      && cur_tetra[3].isXAbove(corners[0], min_height);
        return cur_tetra;
    }

    bool contains(Vec3i64 point) const noexcept
    {
        if (util::vector_arithmetic::any(periodic)) {
            for (int d = 0; d < 3; d++) {
                if (periodic[d] && point[d] < min_dim[d]) {
                    point[d] += box_size[d];
                }
            }
        }
        // Iterate over all tetrahedrons of the domain
        for (int tetra = 0; tetra < 6; tetra++) {
            // Points which are exactly on this plane of the tetrahedron
            // are not accepted to avoid that they are assigned to two domains.
            bool tetraContainsP = tetrahedron_planes[tetra][3].isAbove(point);
            // Iterate over all other sides
            for (int plane = 0; plane < 3; plane++) {
                if (!tetraContainsP) {
                    break;
                }
                tetraContainsP
                    = tetrahedron_planes[tetra][plane].isAboveOrEqual(point);
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

Octagon::Octagon(const std::array<Vec3d, 8> &vertices, double max_cutoff)
    : oi(
        std::make_unique<_Octagon_Impl>(integerizedArray(vertices), max_cutoff))
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
