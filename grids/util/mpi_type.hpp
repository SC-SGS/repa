
#pragma once

#include <cstdint>
#include <mpi.h>

namespace repa {
namespace util {

template <typename T>
struct mpi_type {};

#define TYPE_ASSOC_LIST \
    MPI_TYPE_ASSOC_FUNC(int8_t, MPI_INT8_T) \
    MPI_TYPE_ASSOC_FUNC(int16_t, MPI_INT16_T) \
    MPI_TYPE_ASSOC_FUNC(int32_t, MPI_INT32_T) \
    MPI_TYPE_ASSOC_FUNC(int64_t, MPI_INT64_T) \
    MPI_TYPE_ASSOC_FUNC(uint8_t, MPI_UINT8_T) \
    MPI_TYPE_ASSOC_FUNC(uint16_t, MPI_UINT16_T) \
    MPI_TYPE_ASSOC_FUNC(uint32_t, MPI_UINT32_T) \
    MPI_TYPE_ASSOC_FUNC(uint64_t, MPI_UINT64_T) \
    MPI_TYPE_ASSOC_FUNC(float, MPI_FLOAT) \
    MPI_TYPE_ASSOC_FUNC(double, MPI_DOUBLE)

#define DECLARE_MPI_TYPE(typ) \
    template <> \
    struct mpi_type<typ> { \
        static const MPI_Datatype type; \
    }

#undef MPI_TYPE_ASSOC_FUNC
#define MPI_TYPE_ASSOC_FUNC(typ, mpi_typ) DECLARE_MPI_TYPE(typ);

// All declarations
TYPE_ASSOC_LIST

#undef MPI_TYPE_ASSOC_FUNC

// Convenience macros using decltype
#define MPI_DECLTYPE_T(x) repa::util::mpi_type<decltype(x)>::type
#define MPI_ELEM_DECLTYPE_T(x) repa::util::mpi_type<decltype(x)::value_type>::type

}
}
