
#pragma once

#include <array>

namespace repa {


template <typename T>
using Vec3 = std::array<T, 3>;

typedef Vec3<int> Vec3i;
typedef Vec3<double> Vec3d;

}