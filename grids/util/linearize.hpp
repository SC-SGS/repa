
#pragma once

#include <array>
#include "common_types.hpp"

namespace repa {
namespace util {

template <typename Ret, typename T1, typename T2>
constexpr Ret linearize(const T1 *c, const T2 *grid)
{
  // Cast in case "Ret" is a type capable of holding larger values than "T1" or "T2".
  return (static_cast<Ret>(c[0]) * grid[1] + c[1]) * grid[2] + c[2];
}


template <typename Ret, typename T1, typename T2>
constexpr Ret linearize(const Vec3<T1>& c, const Vec3<T2>& grid)
{
  return linearize<Ret>(c.data(), grid.data());
}

template <typename T>
constexpr T linearize(const T *c, const T *grid)
{
  return ((c[0]) * grid[1] + c[1]) * grid[2] + c[2];
}

template <typename T>
constexpr T linearize(const Vec3<T>& c, const Vec3<T>& grid)
{
  return linearize<T>(c.data(), grid.data());
}

template <typename T>
constexpr Vec3<T> unlinearize(T cidx, const Vec3<T>& grid)
{
  return {{ (cidx / grid[2]) / grid[1],
            (cidx / grid[2]) % grid[1],
            cidx % grid[2] }};
}

}
}

