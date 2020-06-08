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

set(repa_SRC
    grids/glomethod.cpp
    grids/diffusion.cpp
    grids/psdiffusion.cpp
    grids/diff_variants.cpp
    grids/cart.cpp
    grids/gridbased.cpp
    grids/util/initial_partitioning.cpp
    grids/util/neighbor_offsets.cpp
    grids/util/tetra.cpp
    grid_types.cpp
    grid_variants.cpp
    pargrid_factory.cpp
    pargrid.cpp
)

if(WITH_KDPART AND KDPART_FOUND)
  set(repa_SRC
      ${repa_SRC}
      grids/kd_tree.cpp)
endif(WITH_KDPART AND KDPART_FOUND)

if(WITH_P4EST AND P4EST_FOUND)
  set(repa_SRC
      ${repa_SRC}
      grids/p4est.cpp)
endif(WITH_P4EST AND P4EST_FOUND)

if(WITH_PARMETIS AND PARMETIS_FOUND)
  set(repa_SRC
      ${repa_SRC}
      grids/graph.cpp
      grids/hybrid-gp-diff.cpp)
endif(WITH_PARMETIS AND PARMETIS_FOUND)

set(repa_HDR
    repa.hpp
    pargrid.hpp
    pargrid_factory.hpp
    common_types.hpp
    grid_types.hpp
    grid_variants.hpp
    grids/util/range.hpp
    grids/util/simple_variant.hpp
    grids/util/span.hpp
    grids/util/strong_alias.hpp
)

add_library(repa SHARED ${repa_SRC})

if(WITH_KDPART AND KDPART_FOUND)
  target_include_directories(repa PRIVATE ${KDPART_INCLUDE_DIR})
  target_link_libraries(repa PRIVATE ${KDPART_LIBRARIES})
  add_definitions(-DHAVE_KDPART)
endif(WITH_KDPART AND KDPART_FOUND)

if(WITH_P4EST AND P4EST_FOUND)
  target_include_directories(repa PRIVATE ${P4EST_INCLUDE_DIR})
  target_link_libraries(repa PRIVATE ${P4EST_LIBRARIES})
  add_definitions(-DHAVE_P4EST)
endif(WITH_P4EST AND P4EST_FOUND)

if(WITH_PARMETIS AND PARMETIS_FOUND)
  target_include_directories(repa
                             PRIVATE ${PARMETIS_INCLUDE_DIR}
                                     ${METIS_INCLUDE_DIR})
  target_link_libraries(repa PRIVATE ${PARMETIS_LIBRARIES})
  add_definitions(-DHAVE_PARMETIS)
endif(WITH_PARMETIS AND PARMETIS_FOUND)

target_link_libraries(repa
                      PUBLIC ${MPI_CXX_LIBRARIES}
                             ${Boost_MPI_LIBRARIES}
                             ${Boost_SERIALIZATION_LIBRARIES})
target_include_directories(repa PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(repa
                       PRIVATE ${REPA_DEFAULT_COMPILE_OPTIONS})

install(TARGETS repa
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
# Keep folder structure, so no install for all repa_HDR
foreach(hdr ${repa_HDR})
    get_filename_component(dir ${hdr} DIRECTORY)
    install(FILES ${hdr}
            DESTINATION ${CMAKE_INSTALL_PREFIX}/include/repa/${dir})
endforeach()


if(WITH_TESTS)
set(repa_unit_test_SRC
  grids/util/unit_tests/register_unit_tests.cpp
    grids/util/unit_tests/test_push_back_unique.cpp
    grids/util/unit_tests/test_vector_coerce.cpp
    grids/util/unit_tests/test_get_keys.cpp
    grids/util/unit_tests/test_ioptional.cpp
    grids/util/unit_tests/test_linearize.cpp
    grids/util/unit_tests/test_neighbor_offsets.cpp
    grids/util/unit_tests/test_set_union.cpp
    grids/util/unit_tests/test_simple_variant.cpp
    grids/util/unit_tests/test_span.cpp
    grids/util/unit_tests/test_strong_alias.cpp
)

add_library(repa_unit_test SHARED ${repa_unit_test_SRC})
target_link_libraries(repa_unit_test
                      PUBLIC repa)
target_compile_options(repa_unit_test
                      PRIVATE ${REPA_DEFAULT_COMPILE_OPTIONS})
endif(WITH_TESTS)
