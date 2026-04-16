# CMake generated Testfile for 
# Source directory: /data/lsl/SpGEMM/kokkos-kernels/ode/unit_test
# Build directory: /data/lsl/SpGEMM/kokkos-kernels/ode/unit_test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(ode_serial "/data/lsl/SpGEMM/kokkos-kernels/ode/unit_test/KokkosKernels_ode_serial")
set_tests_properties(ode_serial PROPERTIES  _BACKTRACE_TRIPLES "/data/lsl/SpGEMM/kokkos-kernels/cmake/fake_tribits.cmake;212;add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;194;kokkoskernels_add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;146;kokkoskernels_add_executable_and_test;/data/lsl/SpGEMM/kokkos-kernels/ode/unit_test/CMakeLists.txt;66;kokkoskernels_add_unit_test;/data/lsl/SpGEMM/kokkos-kernels/ode/unit_test/CMakeLists.txt;0;")
add_test(ode_openmp "/data/lsl/SpGEMM/kokkos-kernels/ode/unit_test/KokkosKernels_ode_openmp")
set_tests_properties(ode_openmp PROPERTIES  _BACKTRACE_TRIPLES "/data/lsl/SpGEMM/kokkos-kernels/cmake/fake_tribits.cmake;212;add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;194;kokkoskernels_add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;146;kokkoskernels_add_executable_and_test;/data/lsl/SpGEMM/kokkos-kernels/ode/unit_test/CMakeLists.txt;72;kokkoskernels_add_unit_test;/data/lsl/SpGEMM/kokkos-kernels/ode/unit_test/CMakeLists.txt;0;")
