# initialize the testing environment 
# sets up a custom target called build_test
function(TestInit) 
  set(TESTS_NAMES "" PARENT_SCOPE)
endfunction()


# Creates a simple unit test
# - Positional inputs
#     SUBDIR      -   subdirectory containing test files
# - Named inputs
#   - Required
#       SOURCES     -   source files in SUBDIR that are required for compiling the test executable
#   - Optional
#       LIBS        -   *internal* libraries needed for this test
#       INCLUDES    -   *internal* directories needed for this test
#       DRIVER      -   if present, run the test binary under this invocation (e.g. "mpirun -np 4")
#       ARGS        -   arguments to pass to the test executable
# - Named outputs
#   - Optional
#       EXE_OUT     -   absolute path to the test binary created

function(add_unit_test  TEST_SUBDIR TEST_NAME)
  set(options        "")
  set(oneValueArgs   NAMESPACE  EXE_OUT)
  set(multiValueArgs SOURCES LIBS  INCLUDES  DRIVER  ARGS  MPI_NP  OMP_NT)
  cmake_parse_arguments(PARSE_ARGV 1 "TEST" "${options}" "${oneValueArgs}" "${multiValueArgs}")

  #if((NOT TEST_OMP_NT) OR (NOT HAVE_OPENMP))
  #  set(TEST_OMP_NT 1)
  #endif()
  #if((NOT TEST_MPI_NP) OR (NOT HAVE_MPI))
  #  set(TEST_MPI_NP 1)
  #endif()

  message(("TESTNAME = ${TEST_NAME}"))
  set(TEST_EXE ${TEST_NAME})

  # Build the test driver
  set(TEST_DIR ${CMAKE_CURRENT_LIST_DIR}/${TEST_SUBDIR})     # full path to main test dir
  set(TEST_SRC ${TEST_DIR})                                  # directory with source files for this test
  set(TEST_BIN ${CMAKE_CURRENT_BINARY_DIR}/${TEST_SUBDIR})   # where to build this testing binary
  list(TRANSFORM TEST_SOURCES PREPEND ${TEST_SUBDIR}/)
  add_executable(${TEST_EXE} EXCLUDE_FROM_ALL ${TEST_SOURCES})
  set_target_properties(${TEST_EXE}
                        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${TEST_BIN}
                        )
  # set the include directory of the test exe
  target_include_directories(${TEST_EXE} PUBLIC ${CMAKE_SOURCE_DIR})

  if(TEST_INCLUDES)
    target_include_directories(${TEST_EXE} PUBLIC ${TEST_INCLUDES})
  endif()

  # include the LIBS
  if(TEST_LIBS)
    target_link_libraries(${TEST_EXE} PUBLIC ${TEST_LIBS})
  endif()
  
  if(TEST_EXE_OUT)
    # Pass the absolute path to the test driver created to the parent scope
    set(${TEST_EXE_OUT} ${TEST_BIN}/${TEST_EXE} PARENT_SCOPE)
  endif()

  # Base name and command
  set(TEST_CMD  ${TEST_DRIVER} ./${TEST_EXE} ${TEST_ARGS})

  add_test(
    NAME ${TEST_NAME}
    WORKING_DIRECTORY ${TEST_BIN}
    COMMAND ${TEST_CMD}  
  )

  # Append the test names to a common parent_scope variable TESTS_NAMES
  set(TESTS_NAMES ${TESTS_NAMES} ${TEST_NAME} PARENT_SCOPE)

endfunction()