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
#include <type_traits>
#include <variant>

namespace repa {
namespace util {

#define SIMPLE_VARIANT_DEFINE_COMPARISON(_op_)                                 \
    bool operator _op_(const T1 &v)                                            \
    {                                                                          \
        assert(is_first() || is_second());                                     \
        return is_first() && _a _op_ v;                                        \
    }                                                                          \
    bool operator _op_(const T2 &v)                                            \
    {                                                                          \
        assert(is_first() || is_second());                                     \
        return is_second() && _b _op_ v;                                       \
    }                                                                          \
    bool operator _op_(const simple_variant &v)                                \
    {                                                                          \
        assert(is_first() || is_second());                                     \
        assert(v.is_first() || v.is_second());                                 \
        if (is_first())                                                        \
            return v.is_first() && _a _op_ v._a;                               \
        else if (is_second())                                                  \
            return v.is_second() && _b _op_ v._b;                              \
        else                                                                   \
            return false;                                                      \
    }

/** Simple variant of two distinct types.
 * Internally, consists of a union and an additional tag.
 */
template <typename T1,
          typename T2,
          typename = std::enable_if_t<!std::is_same_v<T1, T2>>>
struct simple_variant {

    /** Creates an uninitialized variant
     */
    simple_variant() : _tag(Tag::UNINITIALIZED)
    {
    }

    simple_variant(const simple_variant &other)
    {
        // See comment in operator= to why we to this
        *this = other;
    }

    simple_variant(simple_variant &&other)
    {
        *this = std::forward<simple_variant>(other);
    }

    simple_variant(const T1 &v) : _a(v), _tag(Tag::TYPE_1)
    {
    }
    simple_variant(T1 &&v) : _a(std::forward<T1>(v)), _tag(Tag::TYPE_1)
    {
    }
    simple_variant(const T2 &v) : _b(v), _tag(Tag::TYPE_2)
    {
    }
    simple_variant(T2 &&v) : _b(std::forward<T2>(v)), _tag(Tag::TYPE_2)
    {
    }

    simple_variant &operator=(const simple_variant &other)
    {
        // Don't just copy the union. T1 and T2 might have non-trivial
        // constructors
        if (other.is_first()) {
            *this = other._a;
        }
        else if (other.is_second()) {
            *this = other._b;
        }
        else {
            assert(false);
        }
        return *this;
    }

    simple_variant &operator=(simple_variant &&other)
    {
        if (other.is_first()) {
            *this = std::forward<T1>(other._a);
        }
        else if (other.is_second()) {
            *this = std::forward<T2>(other._b);
        }
        else {
            assert(false);
        }
        return *this;
    }

    explicit operator T1() const
    {
        assert(is_first());
        return _a;
    }
    explicit operator T2() const
    {
        assert(is_second());
        return _b;
    }

    simple_variant &operator=(const T1 &v)
    {
        _tag = Tag::TYPE_1;
        _a = v;
        return *this;
    }

    simple_variant &operator=(T1 &&v)
    {
        _tag = Tag::TYPE_1;
        _a = std::forward<T1>(v);
        return *this;
    }

    simple_variant &operator=(const T2 &v)
    {
        _tag = Tag::TYPE_2;
        _b = v;
        return *this;
    }

    simple_variant &operator=(T2 &&v)
    {
        _tag = Tag::TYPE_2;
        _b = std::forward<T2>(v);
        return *this;
    }

    SIMPLE_VARIANT_DEFINE_COMPARISON(==)
    SIMPLE_VARIANT_DEFINE_COMPARISON(!=)
    SIMPLE_VARIANT_DEFINE_COMPARISON(<)
    SIMPLE_VARIANT_DEFINE_COMPARISON(<=)
    SIMPLE_VARIANT_DEFINE_COMPARISON(>)
    SIMPLE_VARIANT_DEFINE_COMPARISON(>=)

    /** Returns true iff the object holds an element of the first template
     * parameter type. Consider using is() to make your code clearer.
     */
    bool is_first() const
    {
        return _tag == Tag::TYPE_1;
    }

    /** Returns true iff the object holds an element of the second template
     * parameter type. Consider using is() to make your code clearer.
     */
    bool is_second() const
    {
        return _tag == Tag::TYPE_2;
    }

    /** Returns true iff the object holds an element of type "CheckType".
     */
    template <typename CheckType>
    bool is() const
    {
        if constexpr (std::is_same_v<T1, CheckType>)
            return is_first();
        else if constexpr (std::is_same_v<T2, CheckType>)
            return is_second();
        else
            return false;
    }

    /** Returns the held object as "GetType".
     * On debug builds, this function ensures that the held type is actually
     * "GetType".
     */
    template <typename GetType>
    GetType as() const
    {
        static_assert(
            std::is_same_v<T1, GetType> || std::is_same_v<T2, GetType>);
        if constexpr (std::is_same_v<T1, GetType>) {
            if (is_first())
                return _a;
            else
                assert(false);
        }
        else if constexpr (std::is_same_v<T2, GetType>) {
            if (is_second())
                return _b;
            else
                assert(false);
        }
        // Not reached. (see static_assert above)
        ensure_not_reached();
    }

    /** Calls "f" with the held object, if it is of type "T".
     */
    template <typename T, typename F>
    void visit_if(F &&f)
    {
        visitor<T, F>{}(*this, std::forward<F>(f));
    }

    /** @see visit_if(F &&)
     */
    template <typename T, typename F>
    void visit_if(F &&f) const
    {
        visitor<T, F>{}(*this, std::forward<F>(f));
    }

    /** Calls "f1" with the held object, if it is of type "T1" and "f2" if it is
     * of type "T2".
     */
    template <typename F1, typename F2>
    void visit(F1 &&f1, F2 &&f2)
    {
        visit_if<T1>(f1);
        visit_if<T2>(f2);
    }

    /** @see visit(F1, F2)
     */
    template <typename F1, typename F2>
    void visit(F1 &&f1, F2 &&f2) const
    {
        visit_if<T1>(f1);
        visit_if<T2>(f2);
    }

    /** Returns f(obj) where obj is the held object.
     * "f" must be callable with "T1" as well as "T2" as argument.
     */
    template <typename F>
    auto fmap(F &&f)
    {
        if (is_first())
            return f(_a);
        else
            return f(_b);
    }

    /** Returns a std::variant holding a copy of this object.
     */
    std::variant<T1, T2> as_std_variant() const
    {
        std::variant<T1, T2> v;
        auto assign_value = [&v](const auto val) { v = val; };
        visit(assign_value, assign_value);
        assert(std::holds_alternative<T1>(v) || std::holds_alternative<T2>(v));
        return v;
    }

private:
    union {
        T1 _a;
        T2 _b;
    };
    enum class Tag : int_fast8_t {
        UNINITIALIZED = INT8_C(-1),
        TYPE_1 = INT8_C(0),
        TYPE_2 = INT8_C(1)
    } _tag;

    template <typename T, typename F>
    struct visitor;

    template <typename F>
    struct visitor<T1, F> {
        void operator()(simple_variant &v, F &&f)
        {
            if (v.is_first())
                f(v._a);
        }

        void operator()(const simple_variant &v, F &&f)
        {
            if (v.is_first())
                f(v._a);
        }
    };

    template <typename F>
    struct visitor<T2, F> {
        void operator()(simple_variant &v, F &&f)
        {
            if (v.is_second())
                f(v._b);
        }

        void operator()(const simple_variant &v, F &&f)
        {
            if (v.is_second())
                f(v._b);
        }
    };
};

} // namespace util
} // namespace repa

/** Hashes a simple_variant.
 */
template <typename T1, typename T2>
struct std::hash<repa::util::simple_variant<T1, T2>> {
    auto operator()(const repa::util::simple_variant<T1, T2> &val) const
    {
        return _hasher(val.as_std_variant());
    }

private:
    /** Simplest possible implementation: Use std::variant's hash function.
     */
    std::hash<std::variant<T1, T2>> _hasher;
};
