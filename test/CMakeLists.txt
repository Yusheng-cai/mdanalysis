include(TestSetup)

set(_core_libs mdaOrderParameters mdaparallel mdatools mdaxdr)

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

add_unit_test("TestTools/TestCmdArguments"
             "TestCmdArguments"
             SOURCES main.cpp 
             LIBS ${_core_libs}
             ARGS -number 1 2 3 -str a b c -flag true true false)

add_unit_test("TestXdr/"
             "TestTop"
             SOURCES main.cpp 
             LIBS ${_core_libs}
             ARGS ${CMAKE_CURRENT_LIST_DIR}/SampleTrajectory/5CB.itp)

# Now add the test which tests the main program (only test if they run) 
set(OMP_THREADS 1 4 8)
set(OMP_TEST "")
list(APPEND OMP_TEST "TestOP/TestIndus/box")
list(APPEND OMP_TEST "TestOP/TestIndus/sphere")
list(APPEND OMP_TEST "TestOP/TestP2")
list(APPEND OMP_TEST "TestOP/TestP2tilde")
list(APPEND OMP_TEST "TestCalculation/TestRdfResidue")
list(APPEND OMP_TEST "TestCalculation/Testgcost")

set(PROGRAM_PATH "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/mdanalysis")
set(TEST_RUN_SCRIPT "${CMAKE_CURRENT_LIST_DIR}/run_test.sh")

# Find bash
find_program(BASH bash)
if(NOT BASH)
  message(FATAL_ERROR "Unable to find bash for testing")
else()
  message("Found BASH, bash path is ${BASH}")
endif()

# Find python
find_package(Python3 COMPONENTS Interpreter Development)
if(Python3_Interpreter_FOUND)
  message("Python3 found, interpreter is ${Python3_EXECUTABLE}")
else()
  message(FATAL_ERROR, "Python not found")
endif()
    

# Add the test and run with the bash script run_test.sh
foreach(TEST_NAME IN LISTS OMP_TEST)
    foreach(THREADS IN LISTS OMP_THREADS)
        set(TEST_SCRIPT "${TEST_NAME}/input.dat")
        set(TEST_DIR "${CMAKE_CURRENT_LIST_DIR}/${TEST_NAME}")
        set(TEST_REF_DIR "${CMAKE_CURRENT_LIST_DIR}/SampleTrajectory")
        set(TEST_CMD ${BASH} ${TEST_RUN_SCRIPT} -p ${PROGRAM_PATH} -d ${TEST_DIR} 
        -o ${CMAKE_CURRENT_BINARY_DIR}/op.out -r ref.out
             -i ${CMAKE_CURRENT_LIST_DIR}/${TEST_SCRIPT} -n ${THREADS} -f ${TEST_REF_DIR})

        # set working directory to be the directory processed currently
        add_test(NAME "${TEST_NAME}::OMP_${THREADS}"
                COMMAND ${TEST_CMD}
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
        ) 
    endforeach()
endforeach()

set(NONOMP_TEST "")
list(APPEND NONOMP_TEST "TestCalculation/TestPcostz")
list(APPEND NONOMP_TEST "TestCalculation/TestPcost")
list(APPEND NONOMP_TEST "TestCalculation/TestDipole")
list(APPEND NONOMP_TEST "TestCalculation/TestFindNear")

foreach(TEST_NAME IN LISTS NONOMP_TEST)
    set(TEST_SCRIPT "${TEST_NAME}/input.dat")
    set(TEST_DIR "${CMAKE_CURRENT_LIST_DIR}/${TEST_NAME}")
    set(TEST_REF_DIR "${CMAKE_CURRENT_LIST_DIR}/SampleTrajectory")
    set(TEST_CMD ${BASH} ${TEST_RUN_SCRIPT} -p ${PROGRAM_PATH} -d ${TEST_DIR} 
    -o ${CMAKE_CURRENT_BINARY_DIR}/op.out -r ref.out
            -i ${CMAKE_CURRENT_LIST_DIR}/${TEST_SCRIPT} -n 1 -f ${TEST_REF_DIR})

    # set working directory to be the directory processed currently
    add_test(NAME "${TEST_NAME}"
            COMMAND ${TEST_CMD}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
    ) 
endforeach()

add_custom_target(build_test COMMAND ${CMKAE_CTEST_COMMAND}
                              DEPENDS ${TESTS_NAMES})
