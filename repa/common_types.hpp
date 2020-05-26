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
#include <cassert>
#include <exception>
#include <functional>
#include <numeric>
#include <sstream>
#include <type_traits>
#include <vector>

#include <boost/serialization/array.hpp>

namespace repa {

#ifndef NDEBUG
#define REPA_ON_DEBUG(stmt) stmt
#else
#define REPA_ON_DEBUG(stmt) ((void)0)
#endif

namespace __ensure_impl {
[[noreturn]] inline void __ensure_fail(const char *expr,
                                       const char *file,
                                       int line,
                                       const char *func,
                                       const char *msg)
{
    std::printf("Unrecoverable error: Condition failed: `%s' in %s:%d "
                "(%s): %s",
                expr, file, line, func, msg);
    std::abort();
}
} // namespace __ensure_impl

/** Assert equivalent that is *always* ensured and not only if NDEBUG is not set
 * (as assert).
 */
#define ensure(expr, msg)                                                      \
    (static_cast<bool>(expr) ? (void)0                                         \
                             : __ensure_impl::__ensure_fail(                   \
                                 #expr, __FILE__, __LINE__, __func__, msg))

namespace _impl_tt {
template <typename T, typename...>
struct head {
    typedef T type;
};
} // namespace _impl_tt

template <typename T,
          typename TypeTag,
          T Min_Val = std::numeric_limits<T>::min(),
          T Max_Val = std::numeric_limits<T>::max(),
          // We currently want to offer this template only for arithmetic types.
          typename = std::enable_if_t<std::is_arithmetic_v<T>>>
struct StrongAlias {
    using value_type = T;
    using type_tag = TypeTag;

    constexpr void assert_value_domain()
    {
        assert(_value <= Max_Val);
        assert(_value >= Min_Val);
    }

    constexpr StrongAlias() : _value(T{0})
    {
    }

    explicit constexpr StrongAlias(T value) : _value(value)
    {
        assert_value_domain();
    }

    template <typename S,
              typename = std::enable_if_t<
                  !std::is_same_v<T, S> && std::is_integral_v<S>>>
    explicit constexpr StrongAlias(S value)
    {
        REPA_ON_DEBUG(typedef std::uint64_t MaxType);
        REPA_ON_DEBUG(typedef std::int64_t MinType);
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
        assert_value_domain();
    }

    constexpr StrongAlias(const StrongAlias &) = default;
    constexpr StrongAlias(StrongAlias &&) = default;
    StrongAlias &operator=(const StrongAlias &) = default;
    StrongAlias &operator=(StrongAlias &&) = default;

    constexpr operator T() const
    {
        return _value;
    }

    // For-loops need this
    StrongAlias operator++()
    {
        ++_value;
        assert_value_domain();
        return *this;
    }

    StrongAlias operator++(int)
    {
        auto tmp = *this;
        _value++;
        assert_value_domain();
        return tmp;
    }

    StrongAlias &operator+=(const StrongAlias &other)
    {
        _value += other._value;
        assert_value_domain();
        return *this;
    }

    StrongAlias &operator-=(const StrongAlias &other)
    {
        _value -= other._value;
        assert_value_domain();
        return *this;
    }

    StrongAlias operator+(const StrongAlias &other) const
    {
        StrongAlias r{*this};
        r += other;
        return r;
    }

    StrongAlias operator-(const StrongAlias &other) const
    {
        StrongAlias r{*this};
        r -= other;
        return r;
    }

    bool operator==(const StrongAlias &other) const
    {
        return _value == other._value;
    }

    bool operator!=(const StrongAlias &other) const
    {
        return _value != other._value;
    }

    bool operator<=(const StrongAlias &other) const
    {
        return _value <= other._value;
    }

    bool operator<(const StrongAlias &other) const
    {
        return _value < other._value;
    }

    bool operator>=(const StrongAlias &other) const
    {
        return _value >= other._value;
    }
    bool operator>(const StrongAlias &other) const
    {
        return _value > other._value;
    }

    friend class boost::serialization::access;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int file_version)
    {
        ar &_value;
    }

private:
    T _value;
};

template <typename T>
struct is_strong_alias_type {
    static const bool value = false;
};

template <typename T, typename Tag>
struct is_strong_alias_type<StrongAlias<T, Tag>> {
    static const bool value = true;
};

/** Base type for Expression Templates in vec_arith.hpp
 */
template <typename T, size_t N, typename Expr>
struct VecExpression {
    constexpr T operator[](size_t i) const
    {
        return static_cast<const Expr &>(*this)[i];
    }
};

/** Behaves like a std::array.
 */
template <typename T, size_t N>
struct Vec : VecExpression<T, N, Vec<T, N>> {
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

    /** Initializes the data to 0.
     */
    constexpr Vec() : m_data({{T(0), T(0), T(0)}})
    {
    }
    constexpr Vec(Vec &&) = default;
    constexpr Vec(const Vec &) = default;
    Vec &operator=(const Vec &) = default;

    template <typename Expr>
    constexpr Vec(const VecExpression<T, N, Expr> &e)
    {
        for (size_type i = 0; i < N; ++i) {
            m_data[i] = e[i];
        }
    }

    constexpr Vec(underlying_type &&arr)
        : m_data(std::forward<underlying_type>(arr))
    {
    }

    constexpr Vec(const underlying_type &arr) : m_data(arr)
    {
    }

    // Only accept this very general template if all arguments are convertible
    // to "T".
    template <typename... Args,
              typename = _impl_tt::head<
                  typename std::enable_if<
                      std::is_convertible<Args, T>::value>::type...,
                  typename std::enable_if<sizeof...(Args) == N>::type>>
    explicit constexpr Vec(Args... values) : m_data({{values...}})
    {
    }

    constexpr Vec(std::initializer_list<T> list)
    {
        assert(list.size() == size());
        std::copy(std::begin(list), std::end(list), begin());
    }

    constexpr const VecExpression<T, N, Vec<T, N>> &as_expr() const
    {
        return *this;
    }

    iterator begin()
    {
        return m_data.begin();
    }
    constexpr const_iterator cbegin() const
    {
        return m_data.cbegin();
    }
    const_iterator begin() const
    {
        return cbegin();
    }
    iterator end()
    {
        return m_data.end();
    }
    constexpr const_iterator cend() const
    {
        return m_data.cend();
    }
    const_iterator end() const
    {
        return cend();
    }

    reverse_iterator rbegin()
    {
        return m_data.rbegin();
    }
    constexpr const_reverse_iterator crbegin() const
    {
        return m_data.crbegin();
    }
    const_reverse_iterator rbegin() const
    {
        return crbegin();
    }
    reverse_iterator rend()
    {
        return m_data.rend();
    }
    constexpr const_reverse_iterator crend() const
    {
        return m_data.crend();
    }
    const_reverse_iterator rend() const
    {
        return crend();
    }

    constexpr bool empty() const
    {
        return m_data.empty();
    }
    constexpr size_type size() const
    {
        return m_data.size();
    }
    constexpr size_type max_size() const
    {
        return m_data.max_size();
    }
    T *data()
    {
        static_assert(N > 0, "Vec::data() requires actual elements");
        return m_data.data();
    }
    const T *data() const
    {
        static_assert(N > 0, "Vec::data() requires actual elements");
        return m_data.data();
    }

    T &operator[](size_type i)
    {
        assert(i >= 0 && i < N);
        return m_data[i];
    }
    const T &operator[](size_type i) const
    {
        assert(i >= 0 && i < N);
        return m_data[i];
    }

    template <typename BinaryOp>
    T foldl(BinaryOp &&f)
    {
        static_assert(N > 0, "Vec::fold() requires actual elements");
        return std::accumulate(std::next(cbegin()), cend(), *cbegin(),
                               std::forward<BinaryOp>(f));
    }

    constexpr bool operator==(const Vec &other) const
    {
        return m_data == other.m_data;
    }
    constexpr bool operator!=(const Vec &other) const
    {
        return m_data != other.m_data;
    }

    friend class boost::serialization::access;

    template <typename Archive>
    void serialize(Archive &ar, const unsigned int)
    {
        ar &m_data;
    }

    const underlying_type &as_array() const
    {
        return m_data;
    }
    underlying_type &as_array()
    {
        return m_data;
    }

private:
    std::array<T, N> m_data;
};

template <typename T>
using Vec3 = Vec<T, 3>;

typedef Vec3<int> Vec3i;
typedef Vec3<double> Vec3d;

typedef std::function<std::vector<double>(void)> CellMetric;
typedef std::function<double(int, int)> CellCellMetric;
typedef std::function<void(void)> Thunk;

/** Type that behaves like an integral POD and can be restricted to a
 * range. Range is tested at construction time. Should not compile to any
 * overhead on reasonable compilers and optimization levels.
 */
template <typename T,
          T min,
          T max,
          typename = typename std::enable_if<std::is_integral<T>::value>::type>
struct IntegralRange {
    typedef T value_type;

    template <typename S,
              typename
              = typename std::enable_if<std::is_integral<S>::value>::type>
    inline IntegralRange(S v) : value(static_cast<T>(v))
    {
#ifndef NDEBUG
        if (!in_bounds(v))
            throw std::domain_error("IntegralRange: Value not in bounds.");
#endif
    }

    template <typename S,
              typename
              = typename std::enable_if<std::is_integral<S>::value>::type>
    inline IntegralRange operator=(S v)
    {
#ifndef NDEBUG
        if (!in_bounds(v))
            throw std::domain_error("IntegralRange: Value not in bounds.");
#endif
        value = static_cast<T>(v);
        return *this;
    }

    inline operator value_type()
    {
        return value;
    }

    template <typename S,
              typename
              = typename std::enable_if<std::is_integral<S>::value>::type>
    static inline bool in_bounds(S v)
    {
        // Evaluate range check on wide base type.
        typedef typename std::common_type<T, S>::type base_type;
        return static_cast<base_type>(v) >= static_cast<base_type>(min)
               && static_cast<base_type>(v) <= static_cast<base_type>(max);
    }

private:
    value_type value;
};

typedef IntegralRange<std::int_fast32_t, 0, 26> fs_neighidx;

} // namespace repa

template <typename T, typename Tag>
struct std::common_type<repa::StrongAlias<T, Tag>, repa::StrongAlias<T, Tag>> {
    typedef repa::StrongAlias<T, Tag> type;
};

template <typename T, typename Tag>
struct std::hash<repa::StrongAlias<T, Tag>> {
    auto operator()(const repa::StrongAlias<T, Tag> &val) const
    {
        return _hasher(static_cast<T>(val));
    }

private:
    std::hash<T> _hasher;
};
