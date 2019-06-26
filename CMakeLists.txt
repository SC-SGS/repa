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
set(CMAKE_CXX_FLAGS "-Wno-cast-function-type -Wno-unused")

option(WITH_KDPART "Build with ORB LB support" ON)
option(WITH_P4EST "Build with SFC LB support" ON)
option(WITH_PARMETIS "Build with graph partitioning LB support" ON)
option(WITH_TETRA "Build with grid-based LB support" ON)
option(WITH_KDPART "Build with kd grid support" ON)

option(WITH_TESTS "Perform unit tests" ON)

if(WITH_KDPART)
  find_package(KDPart REQUIRED)
endif(WITH_KDPART)

if(WITH_TETRA)
  find_package(Tetra REQUIRED)
endif(WITH_TETRA)

if(WITH_P4EST)
  find_package(P4est REQUIRED)
endif(WITH_P4EST)

if(WITH_PARMETIS)
  find_package(ParMETIS REQUIRED)
endif(WITH_PARMETIS)

find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})

list(APPEND REQ_BOOST_PACKAGES mpi serialization)
if(WITH_TESTS)
  list(APPEND REQ_BOOST_PACKAGES unit_test_framework)
endif(WITH_TESTS)
find_package(Boost "1.59.0" REQUIRED ${REQ_BOOST_PACKAGES})

add_subdirectory(repa)

if(WITH_TESTS)
  set(TEST_MAX_NPROC 4
      CACHE STRING "Maximum number of processes to run the tests with")
  enable_testing()
  add_subdirectory(tests)
endif(WITH_TESTS)

include(FeatureSummary)
feature_summary(WHAT ALL)
