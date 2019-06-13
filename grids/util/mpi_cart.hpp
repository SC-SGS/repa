
#pragma once

#include <mpi.h>
#include "common_types.hpp"

namespace repa {
namespace util {

inline Vec3i mpi_cart_get_dims(MPI_Comm comm)
{
  Vec3i dims, _dummy, _dummy2;
  MPI_Cart_get(comm, 3, dims.data(), _dummy.data(), _dummy2.data());
  return dims;
}

inline Vec3i mpi_cart_get_coords(MPI_Comm comm)
{
  Vec3i coords, _dummy, _dummy2;
  MPI_Cart_get(comm, 3, _dummy.data(), _dummy2.data(), coords.data());
  return coords;
}

}
}


