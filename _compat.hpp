#pragma once

// Compatibility with old code.
#define ROUND_ERROR_PREC 1e-14
#include <cstdio>
#include <mpi.h>
#define errexit()                                                              \
    do {                                                                       \
        std::fprintf(stderr, "Invoking errexit() at %s:%i", __FILE__,          \
                     __LINE__);                                                \
        MPI_Abort(MPI_COMM_WORLD, 1);                                          \
    } while (0)

// Currently only full-periodic grids are supported
#define PERIODIC(x) (1)