
set(P4EST_DIR "" CACHE PATH "P4est directory")

find_path(P4EST_INCLUDE_DIR
          NAMES p4est.h sc.h
          HINTS ${P4EST_DIR}
          ENV C_INCLUDE_PATH
          PATH_SUFFIXES include)

find_library(P4EST_LIBRARIES
             NAMES p4est
             HINTS ${P4EST_DIR}
             ENV LIBRARY_PATH
             PATH_SUFFIXES lib)

find_library(SC_LIBRARIES
             NAMES sc
             HINTS ${P4EST_DIR}
             ENV LIBRARY_PATH
             PATH_SUFFIXES lib)

list(APPEND P4EST_LIBRARIES ${SC_LIBRARIES})


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(P4EST
                                  DEFAULT_MSG
                                  P4EST_LIBRARIES
                                  P4EST_INCLUDE_DIR)
