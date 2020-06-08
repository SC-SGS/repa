
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>

/** Force load the shared library before main.
 */
extern int register_unit_tests();
static int _force_registration = register_unit_tests();
