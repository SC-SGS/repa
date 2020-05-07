
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
#include <algorithm>
#include <cmath>
#include <iostream>

namespace repa {
namespace util {
namespace tetra {

using Vec3i64 = Vec3<int64_t>;
using Vertices = std::array<Vec3i64, 8>;

int16_t precision = 10;
Vec3i64 box_size{0, 0, 0};

// Anonymous namespace for internal linkage
namespace {

using namespace vector_arithmetic;

Vec3i64 integerize(const Vec3d &v)
{
    return static_cast_vec<Vec3i64>(v * static_cast<double>(precision));
}

Vertices integerizedArray(const std::array<Vec3d, 8> &vertices)
{
    Vertices intVert;
    for (int i = 0; i < 8; i++) {
        intVert[i] = integerize(vertices[i]);
    }
    return intVert;
}

/**
 * Returns the minimum and maximum value for each dimension.
 * In the first Vector the minimum values are collected, in the second the
 * maximum values.
 */
std::pair<Vec3i64, Vec3i64> min_max_per_dim(const Vertices &vertices)
{
    Vec3i64 min = box_size;
    Vec3i64 max{0, 0, 0};
    for (int d = 0; d < 3; d++) {
        for (Vec3i64 vertex : vertices) {
            min[d] = std::min(vertex[d], min[d]);
            max[d] = std::max(vertex[d], max[d]);
        }
    }
    return std::make_pair(min, max);
}

/**
 * Returns the indexes for the lower and upper points in the
 * requested dimension.
 */
bool check_orientation_for_vertices(Vertices &vertices)
{
    std::vector<int64_t> lower, upper;
    bool valid = true;
    for (int d = 0; d < 3; d++) {
        // Because the mapping is inverse, an inverse bit is choosen.
        // (eg we map vertices[001] to position 011 in grid)
        // Inverting works as follows 001 -> 100 -> 011
        int d_bit = std::pow(2, 2 - d);
        for (int i = 0; i < 8; i++) {
            if ((i & d_bit) == 0)
                upper.push_back(vertices[i][d]);
            else
                lower.push_back(vertices[i][d]);
        }
        int smallest_upper = *std::min_element(upper.begin(), upper.end());
        int greatest_lower = *std::max_element(lower.begin(), lower.end());
        valid = valid && (smallest_upper > greatest_lower);
    }
    return valid;
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
planes_of_tetrahedron(const std::array<Vec3i64, 4> &vertices)
{
    return {Plane({vertices[0], vertices[1], vertices[2]}),
            Plane({vertices[0], vertices[2], vertices[3]}),
            Plane({vertices[0], vertices[3], vertices[1]}),
            Plane({vertices[1], vertices[3], vertices[2]})};
}

} // namespace

void init_tetra(double min_cell_size, Vec3d box_s)
{
    precision = static_cast<int16_t>(10. / min_cell_size);
    box_size
        = util::vector_arithmetic::static_cast_vec<Vec3i64>(integerize(box_s));
}

struct _Octagon_Impl {
private:
    static const std::array<int, 6> vertexOrder;
    std::array<std::array<Plane, 4>, 6> tetrahedron_planes;
    double min_height;
    bool isValid;
    Vec3<bool> periodic;
    Vec3i64 min_dim;

public:
    _Octagon_Impl() = delete;

    _Octagon_Impl(const Vertices &const_vertices, double max_cutoff)
        : min_height(2. * std::sqrt(3.) * max_cutoff),
          isValid(true),
          periodic({false, false, false})
    {
        Vertices vertices = const_vertices;
        shift_vertices_over_boundaries(vertices);
        Vec3i64 start = vertices[0];
        Vec3i64 end = vertices[7];
        Vec3i64 last = vertices[5];
        for (int i = 0; i < 6; i++) {
            Vec3i64 next = vertices[vertexOrder[i]];
            tetrahedron_planes[i] = generate_tetra({start, end, next, last});
            last = next;
        }
    }

    void shift_vertices_over_boundaries(Vertices &vertices)
    {
        using namespace util::vector_arithmetic;
        Vec3i64 min, max;
        std::tie(min, max) = min_max_per_dim(vertices);
        Vec3<bool> shifted_above = max >= box_size;
        Vec3<bool> shifted_below = min < Vec3i64{0, 0, 0};
        if (any(shifted_above && shifted_below)) {
            std::cerr << "Subdomain too large! This subdomain is larger than "
                      << "the domain itself in at least one dimension!";
        }
        if (any(shifted_below)) {
            for (Vec3i64 &vertex : vertices) {
                vertex += static_cast_vec<Vec3i64>(shifted_below) * box_size;
            }
            std::tie(min, max) = min_max_per_dim(vertices);
        }
        periodic = shifted_above || shifted_below;
        min_dim = min;
        isValid = isValid && check_orientation_for_vertices(vertices);
    }

    /** Returns 4 planes that represent the faces of a tetrahedron.
     * Additionally, updates isValid.
     */
    inline std::array<Plane, 4>
    generate_tetra(const std::array<Vec3i64, 4> &vertices)
    {
        auto cur_tetra = planes_of_tetrahedron(vertices);
        if (has_validity_check())
            isValid = isValid && cur_tetra[0].isXAbove(vertices[3], min_height)
                      && cur_tetra[1].isXAbove(vertices[1], min_height)
                      && cur_tetra[2].isXAbove(vertices[2], min_height)
                      && cur_tetra[3].isXAbove(vertices[0], min_height);
        return cur_tetra;
    }

    bool contains(Vec3i64 point) const noexcept
    {
        for (int d = 0; d < 3; d++) {
            if (periodic[d] && point[d] < min_dim[d]) {
                point[d] += box_size[d];
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

const std::array<int, 6> _Octagon_Impl::vertexOrder = {{1, 3, 2, 6, 4, 5}};

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
