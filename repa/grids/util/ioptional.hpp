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

#include "pargrid.hpp"

namespace repa {
namespace util {

namespace impl {

/**
 * Default emptyness values for use with ioptional
 */
template <typename T>
struct default_empty_value;

template <>
struct default_empty_value<rank_type> {
    // Default empty value for ranks: -1
    // Be careful. MPI might assign a meaning to "-1". This value cannot be
    // stored in an ioptional.
    static const rank_type value = rank_type{-1};
};

template <typename T>
inline constexpr T default_empty_value_v = default_empty_value<T>::value;

} // namespace impl

/** Analog to std::optional with in-band signalling of emptyness.
 * Note that the template parameter "EmptyValue" cannot be stored in an
 * ioptional.
 * The empty value used by default is given by impl::default_empty_value.
 *
 * Note that this class is only a drop-in replacement for std::optional
 * if you do not rely on exceptions. If you do rely on optional throwing
 * exceptions on erroneous value() calls, don't use this class.
 */
template <typename T, T EmptyValue = impl::default_empty_value_v<T>>
struct ioptional {

    constexpr ioptional() : _value(EmptyValue)
    {
    }
    constexpr ioptional(const T &t) : _value(t)
    {
        // Guard user from storing value with special meaning
        assert(_value != EmptyValue);
    }
    constexpr ioptional(T &&t) : _value(std::forward<T>(t))
    {
        // Guard user from storing value with special meaning
        assert(_value != EmptyValue);
    }
    constexpr ioptional(const ioptional &) = default;
    constexpr ioptional(ioptional &&) = default;
    constexpr ioptional &operator=(const ioptional &) = default;
    constexpr ioptional &operator=(ioptional &&) = default;

    /** Returns if a value is held.
     */
    constexpr bool has_value() const
    {
        return _value != EmptyValue;
    }

    constexpr operator bool() const
    {
        return has_value();
    }

    /** Checked access only on debug builds.
     */
    constexpr T value() const
    {
        assert(has_value());
        return _value;
    }

    /** @see value()
     */
    constexpr T operator*() const
    {
        return value();
    }

    /** Aborts if no value is held.
     * (In contrast to value(), always checks, also for non-debug builds).
     */
    constexpr T ensure_value() const
    {
        if (!has_value())
            std::abort();
        return _value;
    }

    /** Returns its value if one is held, else "t".
     */
    constexpr T value_or(const T &t) const
    {
        if (has_value())
            return _value;
        else
            return t;
    }

    /** @see above
     */
    constexpr T value_or(T &&t) const
    {
        if (has_value())
            return _value;
        else
            return std::forward<T>(t);
    }

private:
    T _value;
};

} // namespace util
} // namespace repa