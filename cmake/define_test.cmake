# Copyright 2017-2019 Steffen Hirschmann
#
# This file is part of Repa.
#
# Repa is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Repa is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Repa.  If not, see <https://www.gnu.org/licenses/>.
#

# Defines a unit test.
# Usage:
# define_test(NAME name SRC src.cpp [LIBRARIES lib1 lib2...])
#
function(define_test)
    set(options "")
    set(oneValueArgs NAME SRC)
    set(multiValueArgs LIBRARIES)
    cmake_parse_arguments(TEST
                          "${options}"
                          "${oneValueArgs}"
                          "${multiValueArgs}"
                          ${ARGN})

    add_executable(${TEST_NAME} ${TEST_SRC})

    target_link_libraries(${TEST_NAME}
                          PRIVATE repa
                                  Boost::unit_test_framework
                                  Boost::mpi
                                  Boost::serialization
                                  ${TEST_LIBRARIES})

    target_include_directories(${TEST_NAME}
                               PRIVATE ${CMAKE_SOURCE_DIR})
    target_compile_options(${TEST_NAME}
                           PRIVATE ${REPA_DEFAULT_COMPILE_OPTIONS})

    foreach(nproc RANGE 1 ${TEST_MAX_NPROC} 1)
        add_test(NAME "${TEST_NAME}-${nproc}"
                 COMMAND ${MPIEXEC} "--oversubscribe"
                         ${MPIEXEC_NUMPROC_FLAG} ${nproc}
                         ${MPIEXEC_PREFLAGS} ${TEST_NAME} ${MPIEXEC_POSTFLAGS})
    endforeach()
endfunction(define_test)
