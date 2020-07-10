# Copyright 2017-2020 Steffen Hirschmann
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

# Defines an example program
# Usage:
# define_example(NAME name SRC src.cpp [LIBRARIES lib1 lib2...])
#
function(define_example)
    set(oneValueArgs NAME SRC)
    set(multiValueArgs LIBRARIES)
    cmake_parse_arguments(EXAMPLE
                          "${options}"
                          "${oneValueArgs}"
                          "${multiValueArgs}"
                          ${ARGN})

    add_executable(${EXAMPLE_NAME} EXCLUDE_FROM_ALL ${EXAMPLE_SRC})
    add_dependencies(examples ${EXAMPLE_NAME})

    target_link_libraries(${EXAMPLE_NAME}
                          PRIVATE repa
                                  Boost::mpi
                                  Boost::serialization)

    target_include_directories(${EXAMPLE_NAME}
                               PRIVATE ${CMAKE_SOURCE_DIR})
    target_compile_options(${EXAMPLE_NAME}
                           PRIVATE ${REPA_DEFAULT_COMPILE_OPTIONS})
endfunction(define_example)
