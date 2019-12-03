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
#include <sstream>
#include <type_traits>
#include <vector>

#include <boost/serialization/array.hpp>

namespace repa {

// Custom assert that throws. We can, then, catch assertion exceptions
// in the testing framework.
namespace tassert {

struct AssertionException : public std::exception {
    AssertionException(const char *expr,
                       const char *file,
                       int line,
                       const char *func)
    {
        std::stringstream sstr;
        sstr << file << ":" << line << ": " << func << ": Assertion `" << expr
             << "' failed.";
        err = sstr.str();
    }

    virtual const char *what() const noexcept override
    {
        return err.c_str();
    }

private:
    std::string err;
};

[[noreturn]] inline void
__t_assert__fail(const char *expr, const char *file, int line, const char *func)
{
    throw AssertionException{expr, file, line, func};
}
} // namespace tassert

#ifdef NDEBUG
#define t_assert(expr) ((void)0)
#else
#define t_assert(expr)                                                         \
    (static_cast<bool>(expr) ? (void)0                                         \
                             : repa::tassert::__t_assert__fail(                \
                                   #expr, __FILE__, __LINE__, __func__))
#endif

/** Base type for Expression Templates in vec_arith.hpp
 */
template <typename T, size_t N, typename Expr>
struct VecExpression {
    constexpr T operator[](size_t i) const
    {
        return static_cast<const Expr &>(*this)[i];
    }
};

namespace _impl_tt {
template <typename T, typename...>
struct head {
    typedef T type;
};
} // namespace _impl_tt

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
                  typename std::enable_if<sizeof...(Args) == N>>>
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
        t_assert(in_bounds(v));
    }

    template <typename S,
              typename
              = typename std::enable_if<std::is_integral<S>::value>::type>
    inline IntegralRange operator=(S v)
    {
        t_assert(in_bounds(v));
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
