
#pragma once

#include <cmath>
#include <type_traits>
#include "common_types.hpp"

namespace repa {
namespace util {

template <typename T, typename = std::enable_if_t<std::is_floating_point<T>::value>>
T norm2(const T *v)
{
  return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

template <typename T, typename = std::enable_if_t<std::is_floating_point<T>::value>>
T norm2(const Vec3<T>& v)
{
  return norm2(v.data());
}

template <typename T, typename = std::enable_if_t<std::is_floating_point<T>::value>>
T dist2(const Vec3<T>& v, const Vec3<T>& w)
{
  Vec3<T> vw;
  for (typename Vec3<T>::size_type d = 0; d < v.size(); ++d)
    vw[d] = v[d] - w[d];
  return norm2(vw);
}

}
}
