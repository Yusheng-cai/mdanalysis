include(TestSetup)

set(_core_libs OrderParameters parallel tools xdr)

set(TESTS_NAMES "")

# AtomSelection tests
set(AtomSelectionTest "")

add_unit_test("TestAtomSelection/TestAtomIndexSelection"
              "TestAtomIndexSelection"
              SOURCES main.cpp
              LIBS ${_core_libs}
              ARGS -f ${CMAKE_CURRENT_LIST_DIR}/TestAtomSelection/TestAtomIndexSelection/input.dat -t ${CMAKE_CURRENT_LIST_DIR}
)


add_unit_test("TestAtomSelection/TestAtomTypeSelection"
              "TestAtomTypeSelection"
              SOURCES main.cpp
              LIBS ${_core_libs}
              ARGS -f ${CMAKE_CURRENT_LIST_DIR}/TestAtomSelection/TestAtomTypeSelection/input.dat -t ${CMAKE_CURRENT_LIST_DIR}
)

add_custom_target(build_test COMMAND ${CMKAE_CTEST_COMMAND}
                              DEPENDS ${TESTS_NAMES})