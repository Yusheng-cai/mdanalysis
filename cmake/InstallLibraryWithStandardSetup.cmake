# Installs the library with the given name (and any specified headers)
# - Puts everything into an export set named "${LIB_NAME}Targets"
#
# INPUT
# - Positional inputs
#     LIB_NAME   -   name of the library (same as in the 'add_library' declaration)
#
# - Named inputs
#   - Required
#       (none)
#   - Optional
#       INCLUDE_SUBDIR   -   name of subdirectory where headers will be installed (default: ${LIB_NAME})
#       INCLUDES         -   list of header files to install under ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${INCLUDE_SUBDIR}
# 
function(InstallLibraryWithStandardSetup  LIB_NAME)
  # Parse inputs
  set(options        "")
  set(oneValueArgs   INCLUDE_SUBDIR)
  set(multiValueArgs INCLUDES)
  cmake_parse_arguments(PARSE_ARGV 1 "LIB_INSTALL" "${options}" "${oneValueArgs}" "${multiValueArgs}")

  if(NOT LIB_INSTALL_INCLUDE_SUBDIR)
    set(LIB_INSTALL_INCLUDE_SUBDIR ${LIB_NAME})
  endif()

  # Standard install locations
  include(GNUInstallDirs)

  # Install library
  install(TARGETS ${LIB_NAME}
    EXPORT ${LIB_NAME}Targets
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    OBJECTS  DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )

  # Install public headers
  install(
    FILES
      ${LIB_INSTALL_INCLUDES}
    DESTINATION
      ${CMAKE_INSTALL_INCLUDEDIR}/${LIB_INSTALL_INCLUDE_SUBDIR}
    )
endfunction()