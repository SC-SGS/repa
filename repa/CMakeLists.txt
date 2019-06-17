
if(WITH_KDPART AND KDPART_FOUND)
  include_directories(${KDPART_INCLUDE_DIR})
  list(APPEND LIBRARIES ${KDPART_LIBRARIES})
  add_definitions(-DHAVE_KDPART)
endif(WITH_KDPART AND KDPART_FOUND)

if(WITH_TETRA AND TETRA_FOUND)
  include_directories(${TETRA_INCLUDE_DIR})
  list(APPEND LIBRARIES ${TETRA_LIBRARIES} "-lCGAL" "-gmp")
  add_definitions(-DHAVE_TETRA)
endif(WITH_TETRA AND TETRA_FOUND)

if(WITH_P4EST AND P4EST_FOUND)
  include_directories(${P4EST_INCLUDE_DIR})
  list(APPEND LIBRARIES ${P4EST_LIBRARIES})
  add_definitions(-DHAVE_P4EST)
endif(WITH_P4EST AND P4EST_FOUND)

if(WITH_PARMETIS AND PARMETIS_FOUND)
  include_directories(${PARMETIS_INCLUDE_DIR})
  list(APPEND LIBRARIES ${PARMETIS_LIBRARIES})
  add_definitions(-DHAVE_METIS)
endif(WITH_PARMETIS AND PARMETIS_FOUND)


set(repa_SRC
    grids/graph.cpp
    grids/gridbased.cpp
    grids/diffusion.cpp
    grids/p4est.cpp
    grids/kd_tree.cpp
    grids/hybrid-gp-diff.cpp
    grids/cart.cpp
    grids/util/mpi_type.cpp
    repart/metric.cpp
    pargrid_factory.cpp
    pargrid.cpp
)

set(repa_HDR
    pargrid.hpp
    pargrid_factory.hpp
    common_types.hpp
    generic_dd_grid_types.hpp
    repart/metric.hpp
)

add_library(repa SHARED ${repa_SRC})
target_link_libraries(repa ${MPI_CXX_LIBRARIES})
target_include_directories(repa PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS repa
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
# Keep folder structure, so no install for all repa_HDR
foreach(hdr ${repa_HDR})
    get_filename_component(dir ${hdr} DIRECTORY)
    install(FILES ${hdr} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/repa/${dir})
endforeach()
