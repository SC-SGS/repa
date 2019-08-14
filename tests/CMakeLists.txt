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

include(define_test)

define_test(NAME test_grid_types SRC test_grid_types.cpp)
define_test(NAME test_repart SRC test_repart.cpp)
define_test(NAME test_extra_params SRC test_extra_params.cpp)
define_test(NAME test_cell_numbers SRC test_cell_numbers.cpp)
define_test(NAME test_neighborhood_symmetry SRC test_neighborhood_symmetry.cpp)
define_test(NAME test_process_neighborhood_symmetry SRC test_process_neighborhood_symmetry.cpp)
define_test(NAME test_all_cells_assigned SRC test_all_cells_assigned.cpp)
define_test(NAME test_all_ghost_cells_valid SRC test_all_ghost_cells_valid.cpp)
define_test(NAME test_ghost_exchange_volume SRC test_ghost_exchange_volume.cpp)
define_test(NAME test_ghost_exchange_volume_ext SRC test_ghost_exchange_volume_ext.cpp)
define_test(NAME test_pos_rank_agreement SRC test_pos_rank_agreement.cpp)
