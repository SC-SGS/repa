
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
find_package(Boost "1.53.0" REQUIRED ${REQ_BOOST_PACKAGES})

add_subdirectory(repa)

if(WITH_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif(WITH_TESTS)

include(FeatureSummary)
feature_summary(WHAT ALL)
