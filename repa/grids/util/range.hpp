
#pragma once

#include <boost/range/adaptors.hpp>
#include <boost/range/irange.hpp>

#include "pargrid.hpp"

namespace repa {
namespace util {

/*
template <typename Int>
auto range(Int val) {
    return boost::irange(val);
}
*/

template <typename Int, typename Tag, Int Min, Int Max>
auto range(StrongAlias<Int, Tag, Min, Max> val)
{
    return boost::irange(static_cast<Int>(val))
           | boost::adaptors::transformed(
               [](Int i) { return StrongAlias<Int, Tag, Min, Max>{i}; });
}

} // namespace util
} // namespace repa