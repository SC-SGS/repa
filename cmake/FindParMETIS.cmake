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

set(PARMETIS_DIR "" CACHE PATH "ParMetis directory")

find_path(PARMETIS_INCLUDE_DIR
          parmetis.h
          HINTS ${PARMETIS_DIR}
          ENV C_INCLUDE_PATH
          PATH_SUFFIXES include)

find_path(METIS_INCLUDE_DIR
          metis.h
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
