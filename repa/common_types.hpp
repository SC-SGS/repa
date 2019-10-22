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

#pragma once

#include <array>
#include <functional>
#include <type_traits>
#include <vector>

#include <boost/serialization/array.hpp>

namespace repa {

template <typename T, size_t N>
struct Vec {
    typedef T value_type;
    typedef T *pointer;
    typedef const T *const_pointer;

    typedef T &reference;
    typedef const T &const_reference;

    typedef std::array<T, N> underlying_type;
    typedef typename underlying_type::size_type size_type;
    typedef typename underlying_type::difference_type difference_type;

    typedef typename underlying_type::iterator iterator;
    typedef typename underlying_type::reverse_iterator reverse_iterator;
    typedef typename underlying_type::const_iterator const_iterator;
    typedef
        typename underlying_type::const_reverse_iterator const_reverse_iterator;

    Vec(const underlying_type &a) : a(a)
    {
    }

    template <typename... Args>
    Vec(Args... values) : a({{values...}})
    {
    }

    constexpr iterator begin()
    {
        return a.begin();
    }
    constexpr const_iterator cbegin() const
    {
        return a.cbegin();
    }
    constexpr iterator end()
    {
        return a.end();
    }
    constexpr const_iterator end() const
    {
        return a.cend();
    }

    constexpr bool empty() const
    {
        return a.empty();
    }
    constexpr size_type size() const
    {
        return a.size();
    }
    constexpr size_type max_size() const
    {
        return a.max_size();
    }
    constexpr T *data()
    {
        return a.data();
    }
    constexpr const T *data() const
    {
        return a.data();
    }

    constexpr T &operator[](size_type i)
    {
        return a[i];
    }
    constexpr const T &operator[](size_type i) const
    {
        return a[i];
    }

    bool operator==(const Vec &other) const
    {
        return a == other.a;
    }
    bool operator!=(const Vec &other) const
    {
        return a != other.a;
    }
    bool operator<(const Vec &other) const
    {
        return a < other.a;
    }
    bool operator<=(const Vec &other) const
    {
        return a <= other.a;
    }
    bool operator>(const Vec &other) const
    {
        return a > other.a;
    }
    bool operator>=(const Vec &other) const
    {
        return a >= other.a;
    }

    template <typename Archive>
    void serialize_to(Archive &ar) const
    {
        ar << a;
    }
    template <typename Archive>
    void deserialize_from(Archive &ar)
    {
        ar >> a;
    }

    const underlying_type &as_array() const
    {
        return a;
    }
    underlying_type &as_array()
    {
        return a;
    }

private:
    std::array<T, N> a;
};

template <typename T>
using Vec3 = Vec<T, 3>;

typedef Vec3<int> Vec3i;
typedef Vec3<double> Vec3d;

typedef std::function<std::vector<double>(void)> CellMetric;
typedef std::function<double(int, int)> CellCellMetric;
typedef std::function<void(void)> Thunk;

/** Type that behaves like an integral POD and can be restricted to a range.
 * Range is tested at construction time.
 * Should not compile to any overhead on reasonable compilers and optimization
 * levels.
 */
template <typename T,
          T min,
          T max,
          typename = typename std::enable_if_t<std::is_integral<T>::value>>
struct IntegralRange {
    typedef T value_type;

    template <typename S,
              typename = typename std::enable_if_t<std::is_integral<S>::value>>
    inline IntegralRange(S v) : value(static_cast<T>(v))
    {
        // Evaluate range check on wide base type.
        typedef std::common_type_t<T, S> base_type;
        assert(static_cast<base_type>(v) >= static_cast<base_type>(min)
               && static_cast<base_type>(v) <= static_cast<base_type>(max));
    }
    inline operator value_type()
    {
        return value;
    }

private:
    value_type value;
};

typedef IntegralRange<std::int_fast32_t, 0, 26> fs_neighidx;

} // namespace repa

namespace boost {
namespace serialization {
template <typename Archive, typename T, size_t N>
void load(Archive &ar,
          repa::Vec<T, N> &v,
          const unsigned int /* file_version */)
{
    v.deserialize_from(ar);
}

template <typename Archive, typename T, size_t N>
void save(Archive &ar,
          const repa::Vec<T, N> &v,
          const unsigned int /* file_version */)
{
    v.serialize_to(ar);
}

template <class Archive, typename T, size_t N>
void serialize(Archive &ar, repa::Vec<T, N> &v, const unsigned int file_version)
{
    split_free(ar, v, file_version);
}
} // namespace serialization
} // namespace boost
