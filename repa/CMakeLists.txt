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

set(repa_SRC
    grids/graph.cpp
    grids/gridbased.cpp
    grids/diffusion.cpp
    grids/p4est.cpp
    grids/kd_tree.cpp
    grids/hybrid-gp-diff.cpp
    grids/cart.cpp
    grids/util/neighbor_offsets.cpp
    grids/util/tetra.cpp
    grid_types.cpp
    pargrid_factory.cpp
    pargrid.cpp
)

set(repa_HDR
    repa.hpp
    pargrid.hpp
    pargrid_factory.hpp
    common_types.hpp
    grid_types.hpp
)

add_library(repa SHARED ${repa_SRC})

if(WITH_KDPART AND KDPART_FOUND)
  target_include_directories(repa PUBLIC ${KDPART_INCLUDE_DIR})
  target_link_libraries(repa PUBLIC ${KDPART_LIBRARIES})
  add_definitions(-DHAVE_KDPART)
endif(WITH_KDPART AND KDPART_FOUND)

if(WITH_CGAL AND CGAL_FOUND)
  target_include_directories(repa PUBLIC ${CGAL_INCLUDE_DIR})
  target_link_libraries(repa PUBLIC CGAL::CGAL)
  add_definitions(-DHAVE_CGAL)
endif(WITH_CGAL AND CGAL_FOUND)

if(WITH_P4EST AND P4EST_FOUND)
  target_include_directories(repa PUBLIC ${P4EST_INCLUDE_DIR})
  target_link_libraries(repa PUBLIC ${P4EST_LIBRARIES})
  add_definitions(-DHAVE_P4EST)
endif(WITH_P4EST AND P4EST_FOUND)

if(WITH_PARMETIS AND PARMETIS_FOUND)
  target_include_directories(repa PUBLIC ${PARMETIS_INCLUDE_DIR})
  target_link_libraries(repa PUBLIC ${PARMETIS_LIBRARIES})
  add_definitions(-DHAVE_METIS)
endif(WITH_PARMETIS AND PARMETIS_FOUND)

target_link_libraries(repa 
                      PUBLIC ${MPI_CXX_LIBRARIES}
                             ${Boost_MPI_LIBRARIES}
                             ${Boost_SERIALIZATION_LIBRARIES})
target_include_directories(repa PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS repa
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
# Keep folder structure, so no install for all repa_HDR
foreach(hdr ${repa_HDR})
    get_filename_component(dir ${hdr} DIRECTORY)
    install(FILES ${hdr} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/repa/${dir})
endforeach()
