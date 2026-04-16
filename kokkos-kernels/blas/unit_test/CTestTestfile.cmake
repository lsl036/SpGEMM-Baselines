# CMake generated Testfile for 
# Source directory: /data/lsl/SpGEMM/kokkos-kernels/blas/unit_test
# Build directory: /data/lsl/SpGEMM/kokkos-kernels/blas/unit_test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(blas_serial "/data/lsl/SpGEMM/kokkos-kernels/blas/unit_test/KokkosKernels_blas_serial")
set_tests_properties(blas_serial PROPERTIES  _BACKTRACE_TRIPLES "/data/lsl/SpGEMM/kokkos-kernels/cmake/fake_tribits.cmake;212;add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;194;kokkoskernels_add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;146;kokkoskernels_add_executable_and_test;/data/lsl/SpGEMM/kokkos-kernels/blas/unit_test/CMakeLists.txt;54;kokkoskernels_add_unit_test;/data/lsl/SpGEMM/kokkos-kernels/blas/unit_test/CMakeLists.txt;0;")
add_test(blas_openmp "/data/lsl/SpGEMM/kokkos-kernels/blas/unit_test/KokkosKernels_blas_openmp")
set_tests_properties(blas_openmp PROPERTIES  _BACKTRACE_TRIPLES "/data/lsl/SpGEMM/kokkos-kernels/cmake/fake_tribits.cmake;212;add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;194;kokkoskernels_add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;146;kokkoskernels_add_executable_and_test;/data/lsl/SpGEMM/kokkos-kernels/blas/unit_test/CMakeLists.txt;61;kokkoskernels_add_unit_test;/data/lsl/SpGEMM/kokkos-kernels/blas/unit_test/CMakeLists.txt;0;")
