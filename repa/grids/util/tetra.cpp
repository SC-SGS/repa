
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

/** Constants for a bitset error value.
 */
enum TetraValidityError {
    E_TETRA_OK = 0,      // All good
    E_TETRA_ROTATED = 1, // Wrong orientation
    E_TETRA_INVALID
    = 3 // Invalid node order, not a convex object, not of minimum required size
};

// Anonymous namespace for internal linkage
namespace {

bool _module_initialized = false;
int16_t precision = 10;
Vec3i64 box_size{0, 0, 0};
Vec3d box_size_double{0., 0., 0.};

using namespace vector_arithmetic;

bool is_precision_small_enough(const Vec3d &v)
{
    // We deal with third powers in the code below (isAbove).
    const Vec3d v3 = v * v * v;
    const auto v3_max = std::max(std::max(v3[0], v3[1]), v3[2]);
    const auto v3_min = std::min(std::min(v3[0], v3[1]), v3[2]);

    static const auto i64_max = std::numeric_limits<int64_t>::max();
    static const auto i64_min = std::numeric_limits<int64_t>::min();

    const double p3 = static_cast<double>(precision) * precision * precision;

    // Multiply by ~10 to make room for all additions.
    assert(i64_max / p3 > v3_max * 10);
    assert(i64_min / p3 < v3_min * 10);

    return true;
}

Vec3i64 integerize(const Vec3d &v)
{
    assert(is_precision_small_enough(v));
    return static_cast_vec<Vec3i64>(v * static_cast<double>(precision));
}

static Vertices vertex_coordinates_of_boundingbox(const BoundingBox &bb)
{
    Vec3i64 shift_dim = {0, 0, 0};
    for (size_t i = 0; i < 8; ++i) {
        const auto &vertex = bb.vertices[i];
        const auto &mirror = bb.mirrors[i];
        for (int d = 0; d < 3; d++) {
            // mirrir[d] == 0 if no mirrors were passed to the BoundingBox
            // constructor.
            if (vertex[d] < -mirror[d] * box_size_double[d])
                shift_dim[d] = 1;
        }
    }

    Vertices intVert;
    for (size_t i = 0; i < 8; i++) {
        const auto &vertex = bb.vertices[i];
        const auto &mirror = bb.mirrors[i];

        const Vec3i64 final_shift
            = (shift_dim + static_cast_vec<Vec3i64>(mirror));
        const Vec3i64 shift = box_size * final_shift;

        for (int d = 0; d < 3; ++d) {
            assert(final_shift[d] >= -1);
            assert(final_shift[d] <= 1);

            // Don't truncate (static_cast) negative values.
            if (vertex[d] >= 0.) {
                intVert[i][d] = static_cast<int64_t>(
                                    vertex[d] * static_cast<double>(precision))
                                + shift[d];
            }
            else {
                intVert[i][d]
                    = shift[d]
                      - static_cast<int64_t>(-vertex[d]
                                             * static_cast<double>(precision));
            }
        }
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
    for (const Vec3i64 &vertex : vertices) {
        for (int d = 0; d < 3; d++) {
            min[d] = std::min(vertex[d], min[d]);
            max[d] = std::max(vertex[d], max[d]);
        }
    }
    return std::make_pair(min, max);
}

enum class WhichFunc { MIN, MAX };

/** Returns the minimum and maximum per dimension of a subset of vertices.
 * The subset is chosen via parameter "which_func". It is called with the
 * current vertex number and dimension and returns if it should be taken
 * into account for the maximum or minimum value.
 */
template <typename Pred>
std::pair<Vec3i64, Vec3i64> min_max_per_dim_of_subset(const Vertices &vertices,
                                                      Pred which_func)
{
    constexpr int64_t max_i64 = std::numeric_limits<int64_t>::max();
    Vec3i64 min = {max_i64, max_i64, max_i64};
    Vec3i64 max{0, 0, 0};
    for (size_t i = 0; i < vertices.size(); ++i) {
        for (int d = 0; d < 3; d++) {
            if (which_func(i, d) == WhichFunc::MIN)
                min[d] = std::min(vertices[i][d], min[d]);
            else
                max[d] = std::max(vertices[i][d], max[d]);
        }
    }
    return std::make_pair(min, max);
}

/**
 * Returns the indices for the lower and upper points in the
 * requested dimension.
 */
bool has_valid_orientation(Vertices &vertices)
{
    Vec3i64 biggest_lower, smallest_upper;

    auto max_of_lower__min_of_upper = [](int vertex_no, int dim) {
        // Upper vertex: bit number "dim" is zero.
        if ((vertex_no & (1 << dim)) == 0)
            return WhichFunc::MIN;
        else
            return WhichFunc::MAX;
    };

    std::tie(smallest_upper, biggest_lower)
        = min_max_per_dim_of_subset(vertices, max_of_lower__min_of_upper);

    return smallest_upper[0] > biggest_lower[0]
           && smallest_upper[1] > biggest_lower[1]
           && smallest_upper[2] > biggest_lower[2];
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
        // The distance d of the point to the plane without scaling is:
        // norm(n) * d = <p, n> - height,                                    (1)
        // where <.,.> denotes the scalar product.
        // To be able to use the scaled vectors, we multiply (1) with
        // precision^3.
        //
        // Cross and dot product are bilinear.
        // - normVector as well as its norm are scaled with precision^2
        // - dot(point, normVector) and heightOfPlane with precision^3
        // Hence, dist needs to be scaled with precision^1.

        // Calculate norm as double vector, otherwise variables in function
        // norm() might overflow for too large values of precision.
        const double n = norm(static_cast_vec<Vec3d>(normVector));
        const double d = ceil(dist * precision * n);

        assert(d < static_cast<double>(std::numeric_limits<int64_t>::max()));
        const int64_t int_dist = static_cast<int64_t>(d);

        assert(dot(static_cast_vec<Vec3d>(point),
                   static_cast_vec<Vec3d>(normVector))
               < static_cast<double>(std::numeric_limits<int64_t>::max()));
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
    double prec = 10. / min_cell_size;

    assert(min_cell_size > 0.0);
    assert(prec <= static_cast<double>(std::numeric_limits<int16_t>::max()));

    precision = static_cast<int16_t>(prec);
    box_size = integerize(box_s);
    box_size_double = box_s;
    _module_initialized = true;
}

void init_tetra()
{
    init_tetra(1.0, {1.0, 1.0, 1.0});
}

int16_t get_precision()
{
    assert(_module_initialized);
    return precision;
}

struct _Octagon_Impl {
private:
    static const std::array<int, 6> vertexOrder;
    std::array<std::array<Plane, 4>, 6> tetrahedron_planes;

    /** Minimum height to ensure.
     */
    double min_height;

    /** Bitset of values from TetraValidityError. 0 if this octagon is valid and
     * can be handled by this implementation. Not zero otherwise.
     */
    int valid_status;

    /** True if the octagon has been periodically shifted in this dimension
     */
    Vec3<bool> periodic;

    /** Minimum coordinate in each dimension
     */
    Vec3i64 min_dim;

public:
    _Octagon_Impl() = delete;

    _Octagon_Impl(Vertices vertices, double max_cutoff)
        : min_height(2. * max_cutoff),
          valid_status(E_TETRA_OK),
          periodic({false, false, false})
    {
        initialize_periodicity_handling(vertices);

        // Periodicity habdling requires octagons not to rotate.
        // Otherwise, the periodic shifts in gridbased.cpp are wrong.
        if (!has_valid_orientation(vertices))
            valid_status |= E_TETRA_ROTATED;

        const Vec3i64 start = vertices[0];
        const Vec3i64 end = vertices[7];
        Vec3i64 last = vertices[5];
        for (int i = 0; i < 6; i++) {
            Vec3i64 next = vertices[vertexOrder[i]];
            tetrahedron_planes[i] = generate_tetra({start, end, next, last});
            last = next;
        }
    }

    /** Recognizes if an octagon expands over the periodic boundary and
     * prepares the object to handle contains() queries from particles
     * in different periodic images.
     *
     * Internally:
     * - Shifts all vertices towards the upper periodic boundary
     * - Sets "periodic" flags
     * - Sets the minimum coordinate of any grid point.
     */
    void initialize_periodicity_handling(Vertices &vertices)
    {
        using namespace util::vector_arithmetic;
        Vec3i64 min, max;
        std::tie(min, max) = min_max_per_dim(vertices);
        periodic = max >= box_size;
        min_dim = min;
    }

    /** Returns 4 planes that represent the faces of a tetrahedron.
     * Additionally, updates isValid.
     */
    inline std::array<Plane, 4>
    generate_tetra(const std::array<Vec3i64, 4> &vertices)
    {
        auto cur_tetra = planes_of_tetrahedron(vertices);
        if (has_validity_check()
            && !(cur_tetra[0].isXAbove(vertices[3], min_height)
                 && cur_tetra[1].isXAbove(vertices[1], min_height)
                 && cur_tetra[2].isXAbove(vertices[2], min_height)
                 && cur_tetra[3].isXAbove(vertices[0], min_height)))
            valid_status |= E_TETRA_INVALID;
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
        return valid_status == E_TETRA_OK;
    }
};

const std::array<int, 6> _Octagon_Impl::vertexOrder = {{1, 3, 2, 6, 4, 5}};

// These are declared here because _Octagon_Impl is an imcomplete type in the
// header.
Octagon::Octagon() = default;
Octagon::~Octagon() = default;
Octagon::Octagon(Octagon &&o) = default;

Octagon::Octagon(const BoundingBox &bb, double max_cutoff)
{
    assert(_module_initialized);
    oi = std::make_unique<_Octagon_Impl>(vertex_coordinates_of_boundingbox(bb),
                                         max_cutoff);
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
