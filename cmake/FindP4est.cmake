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
find_package_handle_standard_args(P4est
                                  DEFAULT_MSG
                                  P4EST_LIBRARIES
                                  P4EST_INCLUDE_DIR)
