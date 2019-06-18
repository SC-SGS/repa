
include_directories(${CMAKE_CURRENT_LIST_DIR}/..)

add_executable(test_grid_types test_grid_types.cpp)
target_link_libraries(test_grid_types repa ${Boost_UNIT_TEST_FRAMEWORK_LIBRARIES})
add_test(NAME grid_types COMMAND test_grid_types)

add_executable(test_neighborhood_symmetry test_neighborhood_symmetry.cpp)
target_link_libraries(test_neighborhood_symmetry repa ${Boost_UNIT_TEST_FRAMEWORK_LIBRARIES})
add_test(NAME neighborhood_symmetry COMMAND test_neighborhood_symmetry)