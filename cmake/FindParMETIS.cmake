
set(PARMETIS_DIR "" CACHE PATH "ParMetis directory")

find_path(PARMETIS_INCLUDE_DIR
          parmetis.h
          HINTS ${PARMETIS_DIR}
          ENV C_INCLUDE_PATH
          PATH_SUFFIXES include)

find_library(PARMETIS_LIBRARIES
             parmetis
             HINTS ${PARMETIS_DIR}
             ENV LIBRARY_PATH
             PATH_SUFFIXES lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARMETIS
                                  DEFAULT_MSG
                                  PARMETIS_LIBRARIES
                                  PARMETIS_INCLUDE_DIR)
