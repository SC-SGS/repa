
set(KDPART_DIR "" CACHE PATH "kdpart directory")

find_path(KDPART_INCLUDE_DIR
          kdpart.h
          HINTS ${KDPART_DIR}
          ENV C_INCLUDE_PATH
          PATH_SUFFIXES include)

find_library(KDPART_LIBRARIES
             kdpart
             HINTS ${KDPART_DIR}
             ENV LIBRARY_PATH
             PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(KDPART
                                  DEFAULT_MSG
                                  KDPART_LIBRARIES
                                  KDPART_INCLUDE_DIR)

