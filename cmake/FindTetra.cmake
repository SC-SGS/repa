
set(TETRA_DIR "" CACHE PATH "Tetra directory")

find_path(TETRA_INCLUDE_DIR
          tetra.hpp
          HINTS ${TETRA_DIR}
          ENV C_INCLUDE_PATH
          PATH_SUFFIXES include)

find_library(TETRA_LIBRARIES
             tetra
             HINTS ${TETRA_DIR}
             ENV LIBRARY_PATH
             PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(TETRA
                                  DEFAULT_MSG
                                  TETRA_LIBRARIES
                                  TETRA_INCLUDE_DIR)
