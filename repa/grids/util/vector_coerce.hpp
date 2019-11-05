/**
 * Copyright 2017-2019 Steffen Hirschmann
 * Copyright 2017-2018 Maximilian Wildbrett
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

#pragma once

#include <algorithm>
#include <iterator>
#include <vector>

namespace repa {
namespace util {

namespace impl {

template <typename From, typename To>
struct VectorCoerce;

/**
 * Holds a reference to a vector.
 */
template <typename From>
struct VectorCoerce<From, From> {
    typedef std::vector<From> vec_type;
    typedef const vec_type &vec_cref;
    VectorCoerce(vec_cref in) : out(in)
    {
    }

    vec_cref cref() const
    {
        return out;
    }

    typename vec_type::size_type size() const
    {
        return out.size();
    }

    typename vec_type::const_iterator begin() const
    {
        return out.begin();
    }
    typename vec_type::const_iterator end() const
    {
        return out.end();
    }

private:
    vec_cref out;
};

/**
 * Copies an input vector and casts the elements to type "To".
 */
template <typename From, typename To>
struct VectorCoerce {
    typedef const std::vector<From> &in_vec_cref;
    typedef std::vector<To> vec_type;
    typedef const vec_type &out_vec_cref;

    VectorCoerce(in_vec_cref in)
    {
        out.reserve(in.size());
        std::transform(std::begin(in), std::end(in), std::back_inserter(out),
                       [](const From &v) { return static_cast<To>(v); });
    }

    out_vec_cref cref() const
    {
        return out;
    }

    size_t size() const
    {
        return out.size();
    }

    typename vec_type::const_iterator begin() const
    {
        return out.begin();
    }
    typename vec_type::const_iterator end() const
    {
        return out.end();
    }

private:
    vec_type out;
};

} // namespace impl

/** Coerces the elements of a vector to a different type.
 * Copies, if the source and target type are not the same.
 * Otherwise avoids copying.
 */
template <typename To, typename From>
impl::VectorCoerce<From, To> coerce_vector_to(const std::vector<From> &v)
{
    return impl::VectorCoerce<From, To>{v};
}

} // namespace util
} // namespace repa
