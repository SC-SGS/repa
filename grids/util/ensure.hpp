
#pragma once

#ifndef NDEBUG

#include <cstdio>
#include <sys/signal.h>

// This macro uses kill(2) as a programatic breakpoint.
// (GDB will break at SIGINT.)
#define ENSURE(cond)                                                           \
    do {                                                                       \
        if (!(cond)) {                                                         \
            fprintf(stderr, "Ensure `%s' in %s:%i (%s) failed.\n", #cond,      \
                    __FILE__, __LINE__, __FUNCTION__);                         \
            kill(0, SIGINT);                                                   \
        }                                                                      \
    } while (0)

#else

#define ENSURE(cond)

#endif