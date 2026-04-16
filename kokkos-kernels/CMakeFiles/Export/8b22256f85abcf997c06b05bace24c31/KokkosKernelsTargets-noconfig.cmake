#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Kokkos::kokkoskernels" for configuration ""
set_property(TARGET Kokkos::kokkoskernels APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(Kokkos::kokkoskernels PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_NOCONFIG "CXX"
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libkokkoskernels.a"
  )

list(APPEND _cmake_import_check_targets Kokkos::kokkoskernels )
list(APPEND _cmake_import_check_files_for_Kokkos::kokkoskernels "${_IMPORT_PREFIX}/lib/libkokkoskernels.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
