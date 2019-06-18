
# Defines a unit test.
# Usage:
# define_test(NAME name SRC src.cpp [LIBRARIES lib1 lib2...])
#
function(define_test)
    set(options "")
    set(oneValueArgs NAME SRC)
    set(multiValueArgs LIBRARIES)
    cmake_parse_arguments(TEST
                          "${options}"
                          "${oneValueArgs}"
                          "${multiValueArgs}"
                          ${ARGN})

    add_executable(${TEST_NAME} ${TEST_SRC})

    target_link_libraries(${TEST_NAME}
                          PRIVATE repa
                                  Boost::unit_test_framework
                                  ${TEST_LIBRARIES})

    target_include_directories(${TEST_NAME}
                               PRIVATE ${CMAKE_SOURCE_DIR})

    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endfunction(define_test)
