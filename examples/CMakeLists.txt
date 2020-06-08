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

include(define_example)
add_custom_target(examples)

define_example(NAME example01 SRC example01.cpp)
define_example(NAME example02 SRC example02.cpp)
define_example(NAME example03 SRC example03.cpp)
define_example(NAME example04 SRC example04.cpp)
define_example(NAME example05 SRC example05.cpp)
define_example(NAME example06 SRC example06.cpp)
define_example(NAME example07 SRC example07.cpp)
define_example(NAME example08 SRC example08.cpp)
define_example(NAME example09-seqmd SRC example09-seqmd.cpp)