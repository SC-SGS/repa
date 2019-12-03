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

#include "common_types.hpp"
#include <type_traits>

namespace repa {
namespace util {
namespace vector_arithmetic {

template <typename T, size_t N>
struct Constant : public VecExpression<T, N, Constant<T, N>> {
    constexpr Constant(const T &v) : _v(v)
    {
    }

    constexpr T operator[](size_t i) const
    {
        (void)i;
        return _v;
    }

private:
    const T _v;
};

template <typename T>
Constant<T, 3> constant_vec3(const T &v)
{
    return Constant<T, 3>{v};
}

#define DEFINE_VEC_OP(op_name, op_literal_name, op_operator)                   \
    template <typename T, size_t N, typename Expr1, typename Expr2>            \
    struct op_name : public VecExpression<T, N, op_name<T, N, Expr1, Expr2>> { \
        constexpr op_name(const Expr1 &e1, const Expr2 &e2) : _e1(e1), _e2(e2) \
        {                                                                      \
        }                                                                      \
                                                                               \
        constexpr T operator[](size_t i) const                                 \
        {                                                                      \
            return _e1[i] op_operator _e2[i];                                  \
        }                                                                      \
                                                                               \
    private:                                                                   \
        const Expr1 &_e1;                                                      \
        const Expr2 &_e2;                                                      \
    };                                                                         \
                                                                               \
    /* Vec OP vec */                                                           \
    template <typename T1, typename T2, size_t N, typename Expr1,              \
              typename Expr2>                                                  \
    constexpr op_name<typename std::common_type<T1, T2>::type, N, Expr1,       \
                      Expr2>                                                   \
    operator op_operator(const VecExpression<T1, N, Expr1> &a,                 \
                         const VecExpression<T2, N, Expr2> &b)                 \
    {                                                                          \
        return op_name<typename std::common_type<T1, T2>::type, N, Expr1,      \
                       Expr2>{*static_cast<const Expr1 *>(&a),                 \
                              *static_cast<const Expr2 *>(&b)};                \
    }                                                                          \
                                                                               \
    template <typename T, size_t N, typename Expr1>                            \
    struct op_literal_name                                                     \
        : public VecExpression<T, N, op_literal_name<T, N, Expr1>> {           \
        constexpr op_literal_name(const Expr1 &e1, const T &val)               \
            : _e1(e1), _val(val)                                               \
        {                                                                      \
        }                                                                      \
                                                                               \
        constexpr T operator[](size_t i) const                                 \
        {                                                                      \
            return _e1[i] op_operator _val;                                    \
        }                                                                      \
                                                                               \
    private:                                                                   \
        const Expr1 &_e1;                                                      \
        T _val;                                                                \
    };                                                                         \
                                                                               \
    /* Vec OP constant */                                                      \
    template <typename T1, typename T2, size_t N, typename Expr1,              \
              typename                                                         \
              = typename std::enable_if<std::is_arithmetic<T2>::value>::type>  \
    constexpr op_literal_name<typename std::common_type<T1, T2>::type, N,      \
                              Expr1>                                           \
    operator op_operator(const VecExpression<T1, N, Expr1> &a, const T2 &b)    \
    {                                                                          \
        return op_literal_name<typename std::common_type<T1, T2>::type, N,     \
                               Expr1>{*static_cast<const Expr1 *>(&a), b};     \
    }

DEFINE_VEC_OP(VecSum, VecSumLiteral, +)
DEFINE_VEC_OP(VecSub, VecSubLiteral, -)
DEFINE_VEC_OP(VecMult, VecMultLiteral, *)
DEFINE_VEC_OP(VecDiv, VecDivLiteral, /)
DEFINE_VEC_OP(VecLAnd, VecLAndLiteral, &&)
DEFINE_VEC_OP(VecLOr, VecLOrLiteral, ||)

// Division with literal first (1.0 / x)
template <typename T, size_t N, typename Expr1>
struct LiteralDivVec : public VecExpression<T, N, LiteralDivVec<T, N, Expr1>> {
    constexpr LiteralDivVec(const T &val, const Expr1 &e1) : _val(val), _e1(e1)
    {
    }

    constexpr T operator[](size_t i) const
    {
        return _val / _e1[i];
    }

private:
    T _val;
    const Expr1 &_e1;
};

/* Vec OP constant */
template <typename T1,
          typename T2,
          size_t N,
          typename Expr1,
          typename
          = typename std::enable_if<std::is_arithmetic<T2>::value>::type>
constexpr LiteralDivVec<typename std::common_type<T1, T2>::type, N, Expr1>
operator/(const T2 &b, const VecExpression<T1, N, Expr1> &a)
{
    return LiteralDivVec<typename std::common_type<T1, T2>::type, N, Expr1>{
        b, *static_cast<const Expr1 *>(&a)};
}

#define DEFINE_ASSIGNMENT_OP(op_operator)                                      \
    template <typename T1, typename T2, size_t N, typename Expr>               \
    constexpr Vec<T1, N> &operator op_operator(                                \
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
    constexpr Vec<T1, N> &operator op_operator(Vec<T1, N> &a, const T2 &b)     \
    {                                                                          \
        for (size_t i = 0; i < N; ++i)                                         \
            a[i] op_operator b;                                                \
        return a;                                                              \
    }

DEFINE_ASSIGNMENT_OP(+=)
DEFINE_ASSIGNMENT_OP(-=)
DEFINE_ASSIGNMENT_OP(*=)
DEFINE_ASSIGNMENT_OP(/=)

/** Casting
 */
template <typename ToType, size_t N, typename Expr>
struct VecCast : public VecExpression<ToType, N, VecCast<ToType, N, Expr>> {
    constexpr VecCast(const Expr &e) : _e(e)
    {
    }

    constexpr ToType operator[](size_t i) const
    {
        return static_cast<ToType>(_e[i]);
    }

private:
    const Expr &_e;
};

template <typename ToType, typename FromType, size_t N, typename Expr>
constexpr VecCast<typename ToType::value_type, N, Expr>
static_cast_vec(const VecExpression<FromType, N, Expr> &e)
{
    return VecCast<typename ToType::value_type, N, Expr>{
        *static_cast<const Expr *>(&e)};
}

/** Negation
 */
template <typename T, size_t N, typename Expr>
struct VecNeg : public VecExpression<T, N, VecNeg<T, N, Expr>> {
    constexpr VecNeg(const Expr &e) : _e(e)
    {
    }

    constexpr T operator[](size_t i) const
    {
        return -_e[i];
    }

private:
    const Expr &_e;
};

template <typename T, size_t N, typename Expr>
constexpr VecNeg<T, N, Expr> operator-(const VecExpression<T, N, Expr> &e)
{
    return VecNeg<T, N, Expr>{*static_cast<const Expr *>(&e)};
}

/** Clamping
 */
template <typename T, size_t N, typename Expr, typename ExprLB, typename ExprUB>
struct VecClamp
    : public VecExpression<T, N, VecClamp<T, N, Expr, ExprLB, ExprUB>> {
    constexpr VecClamp(const Expr &e,
                       const ExprLB &lower_bound,
                       const ExprUB &upper_bound)
        : _e(e), _lower_bound(lower_bound), _upper_bound(upper_bound)
    {
    }

    constexpr T operator[](size_t i) const
    {
        return std::min(std::max(_e[i], _lower_bound[i]), _upper_bound[i]);
    }

private:
    const Expr &_e;
    const ExprLB &_lower_bound;
    const ExprUB &_upper_bound;
};

template <typename T, size_t N, typename Expr1, typename Expr2, typename Expr3>
constexpr VecClamp<T, N, Expr1, Expr2, Expr3>
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
    constexpr VecWrap(const Expr &e, const ExprUB &upper_bound)
        : _e(e), _upper_bound(upper_bound)
    {
    }

    constexpr T operator[](size_t i) const
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

template <typename T, size_t N, typename Expr1, typename Expr3>
constexpr VecWrap<T, N, Expr1, Expr3>
vec_wrap(const VecExpression<T, N, Expr1> &v,
         const VecExpression<T, N, Expr3> &upper_bound)
{
    return VecWrap<T, N, Expr1, Expr3>{
        *static_cast<const Expr1 *>(&v),
        *static_cast<const Expr3 *>(&upper_bound)};
}

/** Implements a wrap operation.
 * BE CAREFUL: Result is defined as always positive! (In contrast to modulo.)
 */
template <typename T, size_t N, typename Expr1, typename Expr3>
constexpr VecWrap<T, N, Expr1, Expr3>
operator%(const VecExpression<T, N, Expr1> &v,
          const VecExpression<T, N, Expr3> &upper_bound)
{
    return vec_wrap(v, upper_bound);
}

// Comparisons

#define DEFINE_VEC_COMP(op_name, op_operator)                                  \
    template <typename T, size_t N, typename Expr>                             \
    struct op_name : public VecExpression<bool, N, op_name<T, N, Expr>> {      \
        constexpr op_name(const Expr &e, const T &val) : _e(e), _val(val)      \
        {                                                                      \
        }                                                                      \
                                                                               \
        constexpr bool operator[](size_t i) const                              \
        {                                                                      \
            return _e[i] op_operator _val;                                     \
        }                                                                      \
                                                                               \
    private:                                                                   \
        const Expr &_e;                                                        \
        const T &_val;                                                         \
    };                                                                         \
                                                                               \
    template <typename T, size_t N, typename Expr>                             \
    constexpr op_name<T, N, Expr> operator op_operator(                        \
        const VecExpression<T, N, Expr> &a, const T &val)                      \
    {                                                                          \
        return op_name<T, N, Expr>{*static_cast<const Expr *>(&a), val};       \
    }

DEFINE_VEC_COMP(VecEqualLiteral, ==)
DEFINE_VEC_COMP(VecGreaterLiteral, >)
DEFINE_VEC_COMP(VecGreaterEqualLiteral, >=)
DEFINE_VEC_COMP(VecLessLiteral, <)
DEFINE_VEC_COMP(VecLessEqualLiteral, <=)

template <size_t N, typename Expr>
bool all(const VecExpression<bool, N, Expr> &v)
{
    for (size_t i = 0; i < N; ++i)
        if (!v[i])
            return false;
    return true;
}

template <size_t N, typename Expr>
bool any(const VecExpression<bool, N, Expr> &v)
{
    for (size_t i = 0; i < N; ++i)
        if (v[i])
            return true;
    return false;
}

#define DEFINE_VEC_COMP2(op_name, op_operator)                                 \
    template <size_t N, typename Expr1, typename Expr2>                        \
    struct op_name : public VecExpression<bool, N, op_name<N, Expr1, Expr2>> { \
        constexpr op_name(const Expr1 &e1, const Expr2 &e2) : _e1(e1), _e2(e2) \
        {                                                                      \
        }                                                                      \
                                                                               \
        constexpr bool operator[](size_t i) const                              \
        {                                                                      \
            return _e1[i] op_operator _e2[i];                                  \
        }                                                                      \
                                                                               \
    private:                                                                   \
        const Expr1 &_e1;                                                      \
        const Expr2 &_e2;                                                      \
    };                                                                         \
                                                                               \
    /* Vec OP vec */                                                           \
    template <typename T, size_t N, typename Expr1, typename Expr2>            \
    constexpr op_name<N, Expr1, Expr2> operator op_operator(                   \
        const VecExpression<T, N, Expr1> &a,                                   \
        const VecExpression<T, N, Expr2> &b)                                   \
    {                                                                          \
        return op_name<N, Expr1, Expr2>{*static_cast<const Expr1 *>(&a),       \
                                        *static_cast<const Expr2 *>(&b)};      \
    }

DEFINE_VEC_COMP2(VecEqualVec, ==)
DEFINE_VEC_COMP2(VecGreaterEqualVec, >=)
DEFINE_VEC_COMP2(VecGreaterVec, >)
DEFINE_VEC_COMP2(VecLessEqualVec, <=)
DEFINE_VEC_COMP2(VecLessVec, <)

} // namespace vector_arithmetic
} // namespace util
} // namespace repa