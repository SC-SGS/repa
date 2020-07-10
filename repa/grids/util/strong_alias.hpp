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

#include <cassert>
#include <functional>
#include <limits>
#include <type_traits>

#include <boost/serialization/serialization.hpp>

#ifndef NDEBUG
#define REPA_ON_DEBUG(stmt) (stmt)
#else
#define REPA_ON_DEBUG(stmt) ((void)0)
#endif

namespace repa {
namespace util {

/** Strong alias type for "T".
 * Use a unique, empty struct as "TypeTag".
 *
 * Uninitialized accesses to objects is detected for debug builds.
 *
 * For NDEBUG builds, this class boils down to a simple "T".
 */
template <typename T,
          typename TypeTag,
          // We currently want to offer this template only for arithmetic types.
          typename = std::enable_if_t<std::is_arithmetic_v<T>>>
struct StrongAlias {
    using value_type = T;
    using type_tag = TypeTag;

    /** Marks object as uninitialized.
     */
    constexpr StrongAlias() : _value(T{0})
    {
        REPA_ON_DEBUG(_initialized = false);
    }

    explicit constexpr StrongAlias(T value) : _value(value)
    {
        REPA_ON_DEBUG(_initialized = true);
    }

    /** Enforces that either "S"'s domain is smaller than "T"'s or performs
     * a runtime check for debug builds if "value" fits into "T".
     */
    template <typename S,
              typename = std::enable_if_t<
                  !std::is_same_v<T, S> && std::is_integral_v<S>>>
    explicit constexpr StrongAlias(S value)
    {
#ifndef NDEBUG
        typedef std::uint64_t MaxType;
        typedef std::int64_t MinType;
#endif
        // Either T has a larger domain than S or runtime check if "value" fits
        // into T.
        assert(static_cast<MaxType>(std::numeric_limits<T>::max())
                   >= static_cast<MaxType>(std::numeric_limits<S>::max())
               || static_cast<MaxType>(std::numeric_limits<T>::max())
                      >= static_cast<MaxType>(value));
        assert(static_cast<MinType>(std::numeric_limits<T>::min())
                   <= static_cast<MinType>(std::numeric_limits<S>::min())
               || static_cast<MinType>(std::numeric_limits<T>::min())
                      <= static_cast<MinType>(value));
        _value = static_cast<T>(value);
        REPA_ON_DEBUG(_initialized = true);
    }

    constexpr StrongAlias(const StrongAlias &other) : _value(other._value)
    {
        assert(other._initialized);
        REPA_ON_DEBUG(_initialized = true);
    }

    constexpr StrongAlias(StrongAlias &&other)
        : _value(std::forward<T>(other._value))
    {
        assert(other._initialized);
        REPA_ON_DEBUG(_initialized = true);
    }

    StrongAlias &operator=(StrongAlias other)
    {
        // Copy swap idiom.
        assert(other._initialized);
        REPA_ON_DEBUG(_initialized = true);
        std::swap(_value, other._value);
        return *this;
    }

    constexpr operator T() const
    {
        assert(_initialized);
        return _value;
    }

    constexpr T value() const
    {
        assert(_initialized);
        return _value;
    }

    // For-loops need this
    StrongAlias operator++()
    {
        assert(_initialized);
        ++_value;
        return *this;
    }

    StrongAlias operator++(int)
    {
        assert(_initialized);
        auto tmp = *this;
        _value++;
        return tmp;
    }

    StrongAlias &operator+=(const StrongAlias &other)
    {
        assert(_initialized);
        _value += other._value;
        return *this;
    }

    StrongAlias &operator-=(const StrongAlias &other)
    {
        assert(_initialized && other._initialized);
        _value -= other._value;
        return *this;
    }

    StrongAlias operator+(const StrongAlias &other) const
    {
        assert(_initialized && other._initialized);
        StrongAlias r{*this};
        r += other;
        return r;
    }

    StrongAlias operator-(const StrongAlias &other) const
    {
        assert(_initialized && other._initialized);
        StrongAlias r{*this};
        r -= other;
        return r;
    }

    bool operator==(const StrongAlias &other) const
    {
        assert(_initialized && other._initialized);
        return _value == other._value;
    }

    bool operator!=(const StrongAlias &other) const
    {
        assert(_initialized && other._initialized);
        return _value != other._value;
    }

    bool operator<=(const StrongAlias &other) const
    {
        assert(_initialized && other._initialized);
        return _value <= other._value;
    }

    bool operator<(const StrongAlias &other) const
    {
        assert(_initialized && other._initialized);
        return _value < other._value;
    }

    bool operator>=(const StrongAlias &other) const
    {
        assert(_initialized && other._initialized);
        return _value >= other._value;
    }
    bool operator>(const StrongAlias &other) const
    {
        assert(_initialized && other._initialized);
        return _value > other._value;
    }

    friend class boost::serialization::access;

    template <class Archive>
    void save(Archive &ar, const unsigned int file_version) const
    {
        assert(_initialized);
        ar &_value;
    }

    template <class Archive>
    void load(Archive &ar, const unsigned int file_version)
    {
        ar &_value;
        REPA_ON_DEBUG(_initialized = true);
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

#ifndef NDEBUG
    bool is_initialized() const
    {
        return _initialized;
    }
#endif

private:
    T _value;
#ifndef NDEBUG
    bool _initialized;
#endif
};

} // namespace util
} // namespace repa

/** Implementation of std::hash for StrongAlias types.
 */
template <typename T, typename Tag>
struct std::hash<repa::util::StrongAlias<T, Tag>> {
    auto operator()(const repa::util::StrongAlias<T, Tag> &val) const
    {
        return _hasher(static_cast<T>(val));
    }

private:
    std::hash<T> _hasher;
};
