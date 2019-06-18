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

#include <cstdint>
#include <mpi.h>

namespace repa {
namespace util {

template <typename T>
struct mpi_type {
};

#define TYPE_ASSOC_LIST                                                        \
    MPI_TYPE_ASSOC_FUNC(int8_t, MPI_INT8_T)                                    \
    MPI_TYPE_ASSOC_FUNC(int16_t, MPI_INT16_T)                                  \
    MPI_TYPE_ASSOC_FUNC(int32_t, MPI_INT32_T)                                  \
    MPI_TYPE_ASSOC_FUNC(int64_t, MPI_INT64_T)                                  \
    MPI_TYPE_ASSOC_FUNC(uint8_t, MPI_UINT8_T)                                  \
    MPI_TYPE_ASSOC_FUNC(uint16_t, MPI_UINT16_T)                                \
    MPI_TYPE_ASSOC_FUNC(uint32_t, MPI_UINT32_T)                                \
    MPI_TYPE_ASSOC_FUNC(uint64_t, MPI_UINT64_T)                                \
    MPI_TYPE_ASSOC_FUNC(float, MPI_FLOAT)                                      \
    MPI_TYPE_ASSOC_FUNC(double, MPI_DOUBLE)

#define DECLARE_MPI_TYPE(typ)                                                  \
    template <>                                                                \
    struct mpi_type<typ> {                                                     \
        static const MPI_Datatype type;                                        \
    }

#undef MPI_TYPE_ASSOC_FUNC
#define MPI_TYPE_ASSOC_FUNC(typ, mpi_typ) DECLARE_MPI_TYPE(typ);

// All declarations
TYPE_ASSOC_LIST

#undef MPI_TYPE_ASSOC_FUNC

// Convenience macros using decltype
#define MPI_DECLTYPE_T(x) repa::util::mpi_type<decltype(x)>::type
#define MPI_ELEM_DECLTYPE_T(x)                                                 \
    repa::util::mpi_type<decltype(x)::value_type>::type

} // namespace util
} // namespace repa
