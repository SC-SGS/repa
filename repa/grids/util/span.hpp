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

#include <cstddef> // std::nullptr
#include <iterator>

namespace repa {
namespace util {

/** Minimal version of C++20 span.
 * Stores a non-owning, sized pointer to a continuous sequence of T.
 * 
 * For simplicity:
 * const_span<T> == span<const T>.
 */
template <typename T>
struct span {
    using value_type = T;
    using pointer = T *;
    using const_pointer = std::add_const_t<T> *; // Necessary for boost::range
    using reference = T &;
    using iterator = pointer;
    using const_iterator = const_pointer;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using size_type = std::size_t;
    using difference_type = std::ptrdiff_t;

    constexpr span() : _data(nullptr), _size(0)
    {
    }

    constexpr span(pointer ptr, size_type size) : _data(ptr), _size(size)
    {
    }

    constexpr pointer data() const
    {
        return _data;
    }

    constexpr size_type size() const
    {
        return _size;
    }

    constexpr bool empty() const
    {
        return size() == 0;
    }

    constexpr iterator begin() const
    {
        return _data;
    }

    constexpr iterator end() const
    {
        return _data + _size;
    }

    constexpr reverse_iterator rbegin() const
    {
        return std::make_reverse_iterator(end());
    }

    constexpr reverse_iterator rend() const
    {
        return std::make_reverse_iterator(begin());
    }

    constexpr reference operator[](size_type i) const
    {
        assert(_data);
        assert(i < _size);
        return _data[i];
    }

private:
    pointer _data;
    size_type _size;
};

template <typename T>
using const_span = span<std::add_const_t<T>>;

template <typename T>
constexpr auto make_span(T *ptr, std::size_t n)
{
    return span<T>{ptr, n};
}

template <typename T>
constexpr auto make_span(const T *ptr, std::size_t n)
{
    return span<const T>{ptr, n};
}

template <typename Container>
constexpr auto make_span(Container &cont)
{
    return make_span(cont.data(), cont.size());
}

template <typename Container>
constexpr auto make_span(const Container &cont)
{
    return make_span(cont.data(), cont.size());
}

template <typename T>
constexpr auto make_const_span(T *ptr, std::size_t n)
{
    return span<std::add_const_t<T>>{ptr, n};
}

template <typename T>
constexpr auto make_const_span(const T *ptr, std::size_t n)
{
    return span<const T>{ptr, n};
}

template <typename Container>
constexpr auto make_const_span(Container &cont)
{
    return make_const_span(cont.data(), cont.size());
}

template <typename Container>
constexpr auto make_const_span(const Container &cont)
{
    return make_const_span(cont.data(), cont.size());
}

} // namespace util
} // namespace repa