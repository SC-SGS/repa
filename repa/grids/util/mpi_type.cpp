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

#include "mpi_type.hpp"

namespace repa {
namespace util {

#define DEFINE_MPI_TYPE(typ, mpi_typ)                                          \
    const MPI_Datatype mpi_type<typ>::type = mpi_typ

#undef MPI_TYPE_ASSOC_FUNC
#define MPI_TYPE_ASSOC_FUNC(typ, mpi_typ) DEFINE_MPI_TYPE(typ, mpi_typ);

// Definitions of all "type" statics
TYPE_ASSOC_LIST

#undef MPI_TYPE_ASSOC_FUNC

} // namespace util
} // namespace repa
