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
#pragma once

#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/irange.hpp>

#include "strong_alias.hpp"

namespace repa {
namespace util {

template <typename Int, typename Tag>
auto range(StrongAlias<Int, Tag> val)
{
    return boost::irange(static_cast<Int>(val))
           | boost::adaptors::transformed(
               [](Int i) { return StrongAlias<Int, Tag>{i}; });
}

/** Iota Iterator
 * Iterates through a range integers starting with "value" (see constructor).
 */
template <typename T>
struct iota_iter
    : public boost::iterator_facade<iota_iter<T>,
                                    T,
                                    boost::random_access_traversal_tag,
                                    /* don't use reference type for return */
                                    T> {
private:
    using base_type = boost::
        iterator_facade<iota_iter<T>, T, boost::random_access_traversal_tag, T>;

public:
    using value_type = typename base_type::value_type;
    using difference_type = typename base_type::difference_type;

    iota_iter() = delete;

    iota_iter(T value) : _value(value)
    {
    }

    iota_iter(const iota_iter &) = default;
    iota_iter(iota_iter &&) = default;

    iota_iter &operator=(const iota_iter &) = default;
    iota_iter &operator=(iota_iter &&) = default;

private:
    friend class boost::iterator_core_access;

    value_type dereference() const
    {
        return _value;
    }

    bool equal(const iota_iter &other) const
    {
        return _value == other._value;
    }

    void increment()
    {
        advance(1);
    }

    void decrement()
    {
        advance(-1);
    }

    void advance(difference_type n)
    {
        _value = T{_value + n};
    }

    difference_type distance_to(const iota_iter &other) const
    {
        return other._value - _value;
    }

    T _value;
};

/** Represents an integer range from 0 to N - 1 where N is the parameter "last".
 */
template <typename T>
struct iota_range {
    typedef T value_type;

    iota_iter<T> begin() const
    {
        return iota_iter<T>(T{0});
    }

    iota_iter<T> end() const
    {
        return iota_iter<T>(_last);
    }

    /** Returns the number of integers represented by the range
     */
    size_t size() const
    {
        return _last;
    }

    bool empty() const
    {
        return _last >= 0;
    }

    /** Access a specific element of the range
     * (Although currently not used in repa itself, iota_rage is exposed via the
     * API. This function might be useful to users.)
     */
    T operator[](size_t i) const
    {
        assert(i < size());
        return *std::next(begin(), i);
    }

    iota_range(T last) : _last(last)
    {
    }

private:
    T _last;
};

} // namespace util
} // namespace repa

namespace boost {

/** Allow iota_range to be used in conjunction with boost::range
 */
template <typename T>
struct range_iterator<repa::util::iota_range<T>> {
    typedef repa::util::iota_iter<T> type;
};

template <typename T>
struct range_iterator<const repa::util::iota_range<T>> {
    typedef repa::util::iota_iter<T> type;
};

} // namespace boost