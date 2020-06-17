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

add_library(testenv STATIC "testenv.cpp")
target_include_directories(testenv PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(testenv PRIVATE repa)

# Unit tests
define_test(NAME utest_get_keys SRC unit_tests/test_get_keys.cpp SINGLEPROC)
define_test(NAME utest_integral_range SRC unit_tests/test_integral_range.cpp SINGLEPROC)
define_test(NAME utest_ioptional SRC unit_tests/test_ioptional.cpp SINGLEPROC)
define_test(NAME utest_linearize SRC unit_tests/test_linearize.cpp SINGLEPROC)
define_test(NAME utest_neighbor_offsets SRC unit_tests/test_neighbor_offsets.cpp SINGLEPROC)
define_test(NAME utest_push_back_unique SRC unit_tests/test_push_back_unique.cpp SINGLEPROC)
define_test(NAME utest_set_union SRC unit_tests/test_set_union.cpp SINGLEPROC)
define_test(NAME utest_simple_variant SRC unit_tests/test_simple_variant.cpp SINGLEPROC)
define_test(NAME utest_span SRC unit_tests/test_span.cpp SINGLEPROC)
define_test(NAME utest_strong_alias SRC unit_tests/test_strong_alias.cpp SINGLEPROC)
define_test(NAME utest_tetra SRC unit_tests/test_tetra.cpp SINGLEPROC)
define_test(NAME utest_vec SRC unit_tests/test_vec.cpp SINGLEPROC)
define_test(NAME utest_vec_arith SRC unit_tests/test_vec_arith.cpp SINGLEPROC)
define_test(NAME utest_vector_coerce SRC unit_tests/test_vector_coerce.cpp SINGLEPROC)

# Integration tests
define_test(NAME itest_unknown_exceptions SRC test_unknown_exceptions.cpp SINGLEPROC)
define_test(NAME itest_coloring SRC test_coloring.cpp)
define_test(NAME itest_small_grid SRC test_small_grid.cpp LIBRARIES testenv SINGLEPROC)
define_test(NAME itest_extra_params SRC test_extra_params.cpp LIBRARIES testenv SINGLEPROC)
define_test(NAME itest_cell_numbers SRC test_cell_numbers.cpp LIBRARIES testenv)
define_test(NAME itest_repart_call SRC test_repart_call.cpp LIBRARIES testenv)
define_test(NAME itest_all_cells_uniquely_assigned SRC test_all_cells_uniquely_assigned.cpp LIBRARIES testenv)
define_test(NAME itest_all_ghost_cells_valid SRC test_all_ghost_cells_valid.cpp LIBRARIES testenv)
define_test(NAME itest_cell_neighborhood_symmetry SRC test_cell_neighborhood_symmetry.cpp LIBRARIES testenv)
define_test(NAME itest_process_neighborhood_symmetry SRC test_process_neighborhood_symmetry.cpp LIBRARIES testenv)
define_test(NAME itest_pos_rank_agreement SRC test_pos_rank_agreement.cpp LIBRARIES testenv)
define_test(NAME itest_hybrid_method_switch SRC test_hybrid_method_switch.cpp LIBRARIES testenv)
define_test(NAME itest_ghost_exchange_volume SRC test_ghost_exchange_volume.cpp LIBRARIES testenv)
