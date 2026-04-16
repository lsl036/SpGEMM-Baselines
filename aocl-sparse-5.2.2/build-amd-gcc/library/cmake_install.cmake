# Install script for directory: /data/lsl/SpGEMM/aocl-sparse-5.2.2/library

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so.5.2.2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so.5.2.2")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so.5.2.2"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so.5.2.2")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib" TYPE SHARED_LIBRARY FILES "/data/lsl/SpGEMM/aocl-sparse-5.2.2/build-amd-gcc/library/libaoclsparse.so.5.2.2")
  if(EXISTS "$ENV{DESTDIR}/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so.5.2.2" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so.5.2.2")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so.5.2.2"
         OLD_RPATH "/data/lsl/Software/amd-libflame/lib/LP64:/data/lsl/Software/amd-blis/lib/LP64:/data/lsl/Software/amd-utils/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so.5.2.2")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so")
    file(RPATH_CHECK
         FILE "$ENV{DESTDIR}/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so"
         RPATH "")
  endif()
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib" TYPE SHARED_LIBRARY FILES "/data/lsl/SpGEMM/aocl-sparse-5.2.2/build-amd-gcc/library/libaoclsparse.so")
  if(EXISTS "$ENV{DESTDIR}/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so")
    file(RPATH_CHANGE
         FILE "$ENV{DESTDIR}/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so"
         OLD_RPATH "/data/lsl/Software/amd-libflame/lib/LP64:/data/lsl/Software/amd-blis/lib/LP64:/data/lsl/Software/amd-utils/lib:"
         NEW_RPATH "")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/usr/bin/strip" "$ENV{DESTDIR}/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/lib/libaoclsparse.so")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/include/aoclsparse_types.h;/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/include/aoclsparse_functions.h;/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/include/aoclsparse_convert.h;/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/include/aoclsparse_analysis.h;/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/include/aoclsparse_auxiliary.h;/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/include/aoclsparse_solvers.h;/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/include/aoclsparse.h;/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/include/aoclsparse.hpp;/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/include/aoclsparse_version.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/include" TYPE FILE FILES
    "/data/lsl/SpGEMM/aocl-sparse-5.2.2/library/include/aoclsparse_types.h"
    "/data/lsl/SpGEMM/aocl-sparse-5.2.2/library/include/aoclsparse_functions.h"
    "/data/lsl/SpGEMM/aocl-sparse-5.2.2/library/include/aoclsparse_convert.h"
    "/data/lsl/SpGEMM/aocl-sparse-5.2.2/library/include/aoclsparse_analysis.h"
    "/data/lsl/SpGEMM/aocl-sparse-5.2.2/library/include/aoclsparse_auxiliary.h"
    "/data/lsl/SpGEMM/aocl-sparse-5.2.2/library/include/aoclsparse_solvers.h"
    "/data/lsl/SpGEMM/aocl-sparse-5.2.2/library/include/aoclsparse.h"
    "/data/lsl/SpGEMM/aocl-sparse-5.2.2/library/include/aoclsparse.hpp"
    "/data/lsl/SpGEMM/aocl-sparse-5.2.2/build-amd-gcc/include/aoclsparse_version.h"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/include/kernel-templates/")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/data/lsl/SpGEMM/aocl-sparse-5.2.2/install-amd/include/kernel-templates" TYPE DIRECTORY FILES "/data/lsl/SpGEMM/aocl-sparse-5.2.2/library/src/include/kernel-templates/" FILES_MATCHING REGEX "/[^/]*\\.hpp$" REGEX "/[^/]*\\.h$")
endif()

