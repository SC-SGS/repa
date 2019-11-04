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

project(Repa)

cmake_minimum_required(VERSION 3.9)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(WITH_KDPART "Build with ORB LB support" ON)
option(WITH_P4EST "Build with SFC LB support" ON)
option(WITH_PARMETIS "Build with graph partitioning LB support" ON)
option(WITH_CGAL "Build with grid-based LB support" ON)
option(WITH_KDPART "Build with kd grid support" ON)

option(WITH_TESTS "Perform unit tests" ON)
option(WITH_COVERAGE "Coverage for tests" ON)

if(WITH_KDPART)
  find_package(KDPart REQUIRED)
endif(WITH_KDPART)

if(WITH_P4EST)
  find_package(P4est REQUIRED)
endif(WITH_P4EST)

if(WITH_PARMETIS)
  find_package(ParMETIS REQUIRED)
endif(WITH_PARMETIS)

if(WITH_CGAL)
  set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)
  find_package(CGAL REQUIRED)
endif(WITH_CGAL)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})

list(APPEND REQ_BOOST_PACKAGES mpi serialization)
if(WITH_TESTS)
  list(APPEND REQ_BOOST_PACKAGES unit_test_framework)
endif(WITH_TESTS)
find_package(Boost "1.67.0" REQUIRED ${REQ_BOOST_PACKAGES})

if (Boost_VERSION VERSION_GREATER_EQUAL "1.69" AND Boost_VERSION VERSION_LESS_EQUAL "1.71")
  message(FATAL_ERROR "Boost 1.69 - 1.71 have a bug and are not supported.")
endif()

# Used for librepa itself as well as all tests
# Enable warnings on gcc and clang. Add flags for other compilers below.
set(REPA_DEFAULT_COMPILE_OPTIONS
    $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:GNU>>:
      "-Wall"
      "-Wextra"
      "-Werror"
      "-pedantic"
      "-pedantic-errors"
      #[[ Exclude common things ]]
      "-Wno-unused-parameter"
      "-Wno-sign-compare"
      #[[ The following flags (without "no-") currently cause warnings
          in the code. However, they are not enabled by Wall or Wextra.
          Still list them here as disabled as a reminder. ]]
      "-Wno-conversion"
      "-Wno-sign-conversion" >)

if(WITH_TESTS AND WITH_COVERAGE)
  if (NOT "${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    message(SEND_ERROR "COVERAGE is only possible with Gnu Compiler Collection.")
  else()
    if (NOT "${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
      message(SEND_ERROR "COVERAGE is only compatible with BUILD_TYPE Debug."
                        " Either disable coverage or set BUILD_TYPE to debug.")
    else()
      list(APPEND REPA_DEFAULT_COMPILE_OPTIONS
                  "-g" "--coverage" "-fprofile-arcs" "-ftest-coverage" "-O0")
      set(CMAKE_SHARED_LINKER_FLAGS "--coverage -fprofile-arcs -ftest-coverage")
      set(CMAKE_EXE_LINKER_FLAGS "--coverage -fprofile-arcs -ftest-coverage")

      add_custom_target(coverage
                        COMMAND ${CMAKE_SOURCE_DIR}/util/coverage.sh
                        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
      message(STATUS "Coverage enabled. Use `coverage' target after build to execute coverage tests.")
    endif()
  endif()
endif(WITH_TESTS AND WITH_COVERAGE)

add_subdirectory(repa)

if(WITH_TESTS)
  set(TEST_MAX_NPROC 4
      CACHE STRING "Maximum number of processes to run the tests with")
  enable_testing()
  add_subdirectory(tests)
endif(WITH_TESTS)

include(FeatureSummary)
feature_summary(WHAT ALL)
