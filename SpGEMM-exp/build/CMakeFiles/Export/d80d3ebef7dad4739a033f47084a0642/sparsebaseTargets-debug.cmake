#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "sparsebase::sparsebase" for configuration "Debug"
set_property(TARGET sparsebase::sparsebase APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(sparsebase::sparsebase PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libsparsebase.a"
  )

list(APPEND _cmake_import_check_targets sparsebase::sparsebase )
list(APPEND _cmake_import_check_files_for_sparsebase::sparsebase "${_IMPORT_PREFIX}/lib/libsparsebase.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
