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