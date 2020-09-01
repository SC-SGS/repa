/**
 * Copyright 2017-2019 Steffen Hirschmann
 *
 * This file is part of Repa.

 * Repa is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.

 * Repa is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with Repa.  If not, see <https://www.gnu.org/licenses/>.
 */

#pragma once

#include "common_types.hpp"
#include <cmath>
#include <type_traits>

/** Expression templates for efficient arithmetics based on "Vec".
 *
 * Do *NOT* include them in any header. Also, only once per translation unit.
 */

namespace repa {
namespace util {
namespace vector_arithmetic {

// Specifiers for all operator or function templates returning
// expression template structs as well as methods of expression template
// structs.
#define __REPA__EX_TMPL_EXPORT constexpr inline

/** Expression template functor returning a constant.
 */
template <typename T, size_t N>
struct Constant : public VecExpression<T, N, Constant<T, N>> {
    __REPA__EX_TMPL_EXPORT Constant(const T &v) : _v(v)
    {
    }

    __REPA__EX_TMPL_EXPORT T operator[](size_t i) const
    {
        (void)i;
        return _v;
    }

private:
    const T _v;
};

/** Returns a constant expression of dimension 3.
 */
template <typename T>
__REPA__EX_TMPL_EXPORT Constant<T, 3> constant_vec3(const T &v)
{
    return Constant<T, 3>{v};
}

// Internal definition of Vector-Vector, Literal-Vector and Vector-Literal
// operators. Defines expression template functors and corresponding operators.
// "op_operator" is the operator to define (e.g. `+')
// "op_return_type" is the return type of the operator (either `T' in case of
// arithmetics or `bool' in case of comparisons) "op_name", "op_literal_name",
// "op_literal_first_name" are names for the espression template classes.
#define __REPA__DEFINE_VEC_OP_IMPL(op_operator, op_return_type, op_name,       \
                                   op_literal_name, op_literal_first_name)     \
    /* Vec OP vec */                                                           \
    template <typename T, size_t N, typename Expr1, typename Expr2>            \
    struct op_name : public VecExpression<op_return_type, N,                   \
                                          op_name<T, N, Expr1, Expr2>> {       \
        __REPA__EX_TMPL_EXPORT op_name(const Expr1 &e1, const Expr2 &e2)       \
            : _e1(e1), _e2(e2)                                                 \
        {                                                                      \
        }                                                                      \
                                                                               \
        __REPA__EX_TMPL_EXPORT op_return_type operator[](size_t i) const       \
        {                                                                      \
            return _e1[i] op_operator _e2[i];                                  \
        }                                                                      \
                                                                               \
    private:                                                                   \
        const Expr1 &_e1;                                                      \
        const Expr2 &_e2;                                                      \
    };                                                                         \
                                                                               \
    template <typename T1, typename T2, size_t N, typename Expr1,              \
              typename Expr2>                                                  \
    __REPA__EX_TMPL_EXPORT                                                     \
        op_name<typename std::common_type<T1, T2>::type, N, Expr1, Expr2>      \
        operator op_operator(const VecExpression<T1, N, Expr1> &a,             \
                             const VecExpression<T2, N, Expr2> &b)             \
    {                                                                          \
        return op_name<typename std::common_type<T1, T2>::type, N, Expr1,      \
                       Expr2>{*static_cast<const Expr1 *>(&a),                 \
                              *static_cast<const Expr2 *>(&b)};                \
    }                                                                          \
                                                                               \
    /* Literal OP vec */                                                       \
    template <typename T, size_t N, typename Expr1>                            \
    struct op_literal_first_name                                               \
        : public VecExpression<op_return_type, N,                              \
                               op_literal_first_name<T, N, Expr1>> {           \
        __REPA__EX_TMPL_EXPORT op_literal_first_name(const T &val,             \
                                                     const Expr1 &e1)          \
            : _val(val), _e1(e1)                                               \
        {                                                                      \
        }                                                                      \
                                                                               \
        __REPA__EX_TMPL_EXPORT op_return_type operator[](size_t i) const       \
        {                                                                      \
            return _val op_operator _e1[i];                                    \
        }                                                                      \
                                                                               \
    private:                                                                   \
        T _val;                                                                \
        const Expr1 &_e1;                                                      \
    };                                                                         \
                                                                               \
    template <typename T1, typename T2, size_t N, typename Expr1,              \
              typename                                                         \
              = typename std::enable_if<std::is_arithmetic<T2>::value>::type>  \
    __REPA__EX_TMPL_EXPORT op_literal_first_name<                              \
        typename std::common_type<T1, T2>::type, N, Expr1>                     \
    operator op_operator(const T2 &a, const VecExpression<T1, N, Expr1> &b)    \
    {                                                                          \
        return op_literal_first_name<typename std::common_type<T1, T2>::type,  \
                                     N, Expr1>{                                \
            a, *static_cast<const Expr1 *>(&b)};                               \
    }                                                                          \
                                                                               \
    /* Vec OP literal */                                                       \
    template <typename T, size_t N, typename Expr1>                            \
    struct op_literal_name                                                     \
        : public VecExpression<op_return_type, N,                              \
                               op_literal_name<T, N, Expr1>> {                 \
        __REPA__EX_TMPL_EXPORT op_literal_name(const Expr1 &e1, const T &val)  \
            : _e1(e1), _val(val)                                               \
        {                                                                      \
        }                                                                      \
                                                                               \
        __REPA__EX_TMPL_EXPORT op_return_type operator[](size_t i) const       \
        {                                                                      \
            return _e1[i] op_operator _val;                                    \
        }                                                                      \
                                                                               \
    private:                                                                   \
        const Expr1 &_e1;                                                      \
        T _val;                                                                \
    };                                                                         \
                                                                               \
    template <typename T1, typename T2, size_t N, typename Expr1,              \
              typename                                                         \
              = typename std::enable_if<std::is_arithmetic<T2>::value>::type>  \
    __REPA__EX_TMPL_EXPORT                                                     \
        op_literal_name<typename std::common_type<T1, T2>::type, N, Expr1>     \
        operator op_operator(const VecExpression<T1, N, Expr1> &a,             \
                             const T2 &b)                                      \
    {                                                                          \
        return op_literal_name<typename std::common_type<T1, T2>::type, N,     \
                               Expr1>{*static_cast<const Expr1 *>(&a), b};     \
    }

// Token concatenation
#define __REPA__CONCAT_IMPL3(a, b, c) a##b##c
#define __REPA__CONCAT3(a, b, c) __REPA__CONCAT_IMPL3(a, b, c)

// For vector operator definition use. Automatically provides the required class
// names
#define __REPA__DEFINE_VEC_OP(op_operator, op_return_type, op_name)            \
    __REPA__DEFINE_VEC_OP_IMPL(op_operator, op_return_type,                    \
                               __REPA__CONCAT3(Vec, op_name, Vec),             \
                               __REPA__CONCAT3(Vec, op_name, Literal),         \
                               __REPA__CONCAT3(Literal, op_name, Vec))

// Define all standard operations
__REPA__DEFINE_VEC_OP(+, T, Plus)
__REPA__DEFINE_VEC_OP(-, T, Minus)
__REPA__DEFINE_VEC_OP(*, T, Times)
__REPA__DEFINE_VEC_OP(/, T, Divide)
__REPA__DEFINE_VEC_OP(&&, T, LogiAnd)
__REPA__DEFINE_VEC_OP(||, T, LogiOr)
__REPA__DEFINE_VEC_OP(>>, T, RightShift)
__REPA__DEFINE_VEC_OP(<<, T, LeftShift)
__REPA__DEFINE_VEC_OP(==, bool, Equal)
__REPA__DEFINE_VEC_OP(!=, bool, Inequal)
__REPA__DEFINE_VEC_OP(>=, bool, GreaterEqual)
__REPA__DEFINE_VEC_OP(>, bool, Greater)
__REPA__DEFINE_VEC_OP(<=, bool, LessEqual)
__REPA__DEFINE_VEC_OP(<, bool, Less)

// Define Vec +=/-=/... assignment operators.
#define __REPA__DEFINE_VEC_ASSIGNMENT_OP(op_operator)                          \
    template <typename T1, typename T2, size_t N, typename Expr>               \
    __REPA__EX_TMPL_EXPORT Vec<T1, N> &operator op_operator(                   \
        Vec<T1, N> &a, const VecExpression<T2, N, Expr> &b)                    \
    {                                                                          \
        for (size_t i = 0; i < N; ++i)                                         \
            a[i] op_operator b[i];                                             \
        return a;                                                              \
    }                                                                          \
                                                                               \
    template <typename T1, size_t N, typename T2,                              \
              typename                                                         \
              = typename std::enable_if<std::is_arithmetic<T2>::value>::type>  \
    __REPA__EX_TMPL_EXPORT Vec<T1, N> &operator op_operator(Vec<T1, N> &a,     \
                                                            const T2 &b)       \
    {                                                                          \
        for (size_t i = 0; i < N; ++i)                                         \
            a[i] op_operator b;                                                \
        return a;                                                              \
    }

__REPA__DEFINE_VEC_ASSIGNMENT_OP(+=)
__REPA__DEFINE_VEC_ASSIGNMENT_OP(-=)
__REPA__DEFINE_VEC_ASSIGNMENT_OP(*=)
__REPA__DEFINE_VEC_ASSIGNMENT_OP(/=)

/** Casting expression template operator
 */
template <typename ToType, size_t N, typename Expr>
struct VecCast : public VecExpression<ToType, N, VecCast<ToType, N, Expr>> {
    __REPA__EX_TMPL_EXPORT VecCast(const Expr &e) : _e(e)
    {
    }

    __REPA__EX_TMPL_EXPORT ToType operator[](size_t i) const
    {
        return static_cast<ToType>(_e[i]);
    }

private:
    const Expr &_e;
};

/** Casts all elements of the VecExpression at evaluation time using a
 * static_cast to type "ToType".
 */
template <typename ToType, typename FromType, size_t N, typename Expr>
__REPA__EX_TMPL_EXPORT VecCast<typename ToType::value_type, N, Expr>
static_cast_vec(const VecExpression<FromType, N, Expr> &e)
{
    return VecCast<typename ToType::value_type, N, Expr>{
        *static_cast<const Expr *>(&e)};
}

/** Negation
 */
template <typename T, size_t N, typename Expr>
struct VecNeg : public VecExpression<T, N, VecNeg<T, N, Expr>> {
    __REPA__EX_TMPL_EXPORT VecNeg(const Expr &e) : _e(e)
    {
    }

    __REPA__EX_TMPL_EXPORT T operator[](size_t i) const
    {
        return -_e[i];
    }

private:
    const Expr &_e;
};

/** Negates all elements of "e" at evaluation time.
 */
template <typename T, size_t N, typename Expr>
__REPA__EX_TMPL_EXPORT VecNeg<T, N, Expr>
operator-(const VecExpression<T, N, Expr> &e)
{
    return VecNeg<T, N, Expr>{*static_cast<const Expr *>(&e)};
}

/** Clamping
 */
template <typename T, size_t N, typename Expr, typename ExprLB, typename ExprUB>
struct VecClamp
    : public VecExpression<T, N, VecClamp<T, N, Expr, ExprLB, ExprUB>> {
    __REPA__EX_TMPL_EXPORT VecClamp(const Expr &e,
                                    const ExprLB &lower_bound,
                                    const ExprUB &upper_bound)
        : _e(e), _lower_bound(lower_bound), _upper_bound(upper_bound)
    {
    }

    __REPA__EX_TMPL_EXPORT T operator[](size_t i) const
    {
        return std::min(std::max(_e[i], _lower_bound[i]), _upper_bound[i]);
    }

private:
    const Expr &_e;
    const ExprLB &_lower_bound;
    const ExprUB &_upper_bound;
};

/** Clampls all elements of "v" at evaluation time between "lower_bound" and
 * "upper_bound". I.e. vec_clamp(a, b, c)[i] == min(max(a[i], b[i]), c[i])
 */
template <typename T, size_t N, typename Expr1, typename Expr2, typename Expr3>
__REPA__EX_TMPL_EXPORT VecClamp<T, N, Expr1, Expr2, Expr3>
vec_clamp(const VecExpression<T, N, Expr1> &v,
          const VecExpression<T, N, Expr2> &lower_bound,
          const VecExpression<T, N, Expr3> &upper_bound)
{
    return VecClamp<T, N, Expr1, Expr2, Expr3>{
        *static_cast<const Expr1 *>(&v),
        *static_cast<const Expr2 *>(&lower_bound),
        *static_cast<const Expr3 *>(&upper_bound)};
}

/** Periodic wrapping
 */
template <typename T, size_t N, typename Expr, typename ExprUB>
struct VecWrap : public VecExpression<T, N, VecWrap<T, N, Expr, ExprUB>> {
    __REPA__EX_TMPL_EXPORT VecWrap(const Expr &e, const ExprUB &upper_bound)
        : _e(e), _upper_bound(upper_bound)
    {
    }

    __REPA__EX_TMPL_EXPORT T operator[](size_t i) const
    {
        T val = _e[i];
        const T ub = _upper_bound[i];
        if (val >= ub) {
            while (val >= ub)
                val -= ub;
        }
        else if (val < 0) {
            while (val < 0)
                val += ub;
        }
        return val;
    }

private:
    const Expr &_e;
    const ExprUB &_upper_bound;
};

/** Wraps elements of "v" periodically at evaluation time.
 * @see VecWrap::operator[]
 */
template <typename T, size_t N, typename Expr1, typename Expr3>
__REPA__EX_TMPL_EXPORT VecWrap<T, N, Expr1, Expr3>
vec_wrap(const VecExpression<T, N, Expr1> &v,
         const VecExpression<T, N, Expr3> &upper_bound)
{
    return VecWrap<T, N, Expr1, Expr3>{
        *static_cast<const Expr1 *>(&v),
        *static_cast<const Expr3 *>(&upper_bound)};
}

/** Alias for vec_wrap()
 * BE CAREFUL: Result is defined as always positive! (In contrast to modulo.)
 * @see vec_wrap()
 */
template <typename T, size_t N, typename Expr1, typename Expr3>
__REPA__EX_TMPL_EXPORT VecWrap<T, N, Expr1, Expr3>
operator%(const VecExpression<T, N, Expr1> &v,
          const VecExpression<T, N, Expr3> &upper_bound)
{
    return vec_wrap(v, upper_bound);
}

/** Returns true if all elements of "v" are true.
 */
template <size_t N, typename Expr>
bool all(const VecExpression<bool, N, Expr> &v)
{
    for (size_t i = 0; i < N; ++i)
        if (!v[i])
            return false;
    return true;
}

/** Returns true if any element of "v" is true.
 */
template <size_t N, typename Expr>
bool any(const VecExpression<bool, N, Expr> &v)
{
    for (size_t i = 0; i < N; ++i)
        if (v[i])
            return true;
    return false;
}

/** Cross product
 */
template <typename T, size_t N, typename Expr1, typename Expr2>
struct VecCross : public VecExpression<T, N, VecCross<T, N, Expr1, Expr2>> {
    __REPA__EX_TMPL_EXPORT VecCross(const Expr1 &u, const Expr2 &v)
        : _u(u), _v(v)
    {
    }

    __REPA__EX_TMPL_EXPORT T operator[](size_t i) const
    {
        return _u[(i + 1) % N] * _v[(i + 2) % N]
               - _u[(i + 2) % N] * _v[(i + 1) % N];
    }

private:
    const Expr1 &_u;
    const Expr2 &_v;
};

/** Vector cross product.
 */
template <typename T1, typename T2, size_t N, typename Expr1, typename Expr2>
__REPA__EX_TMPL_EXPORT
    VecCross<typename std::common_type<T1, T2>::type, N, Expr1, Expr2>
    cross(const VecExpression<T1, N, Expr1> &a,
          const VecExpression<T2, N, Expr2> &b)
{
    return VecCross<typename std::common_type<T1, T2>::type, N, Expr1, Expr2>{
        *static_cast<const Expr1 *>(&a), *static_cast<const Expr2 *>(&b)};
}

/** Dot product.
 */
template <typename T1, typename T2, size_t N, typename Expr1, typename Expr2>
__REPA__EX_TMPL_EXPORT typename std::common_type<T1, T2>::type
dot(const VecExpression<T1, N, Expr1> &v1,
    const VecExpression<T2, N, Expr2> &v2)
{
    typename std::common_type<T1, T2>::type result{0};
    for (size_t i = 0; i < N; ++i)
        result += v1[i] * v2[i];
    return result;
}

template <typename T1, size_t N, typename Expr1>
__REPA__EX_TMPL_EXPORT T1 norm2(const VecExpression<T1, N, Expr1> &v)
{
    T1 result{0};
    for (size_t i = 0; i < N; ++i) {
        const T1 val = v[i]; // Evaluate only once
        result += val * val;
    }
    return result;
}

template <typename T1, size_t N, typename Expr1>
__REPA__EX_TMPL_EXPORT double norm(const VecExpression<T1, N, Expr1> &v)
{
    return std::sqrt(norm2(v));
}

/** Sum of the absoulte values of a Vector.
 *  Can be used to compute the L1-norm of the vector.
 */
template <typename T1, size_t N, typename Expr1>
T1 sumOfAbs(const VecExpression<T1, N, Expr1> &v)
{
    T1 result{0};
    for (size_t i = 0; i < N; ++i)
        result += std::abs(v[i]);
    return result;
}

/** Product of the values of a Vector.
 */
template <typename T1, size_t N, typename Expr1>
T1 product(const VecExpression<T1, N, Expr1> &v)
{
    T1 result{1};
    for (size_t i = 0; i < N; ++i)
        result *= v[i];
    return result;
}

} // namespace vector_arithmetic
} // namespace util
} // namespace repa
