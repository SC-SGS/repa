
#include "mpi_type.hpp"

namespace repa {
namespace util {

#define DEFINE_MPI_TYPE(typ, mpi_typ) \
    const MPI_Datatype mpi_type<typ>::type = mpi_typ

#undef MPI_TYPE_ASSOC_FUNC
#define MPI_TYPE_ASSOC_FUNC(typ, mpi_typ) DEFINE_MPI_TYPE(typ, mpi_typ);

// Definitions of all "type" statics
TYPE_ASSOC_LIST

#undef MPI_TYPE_ASSOC_FUNC

}
}
