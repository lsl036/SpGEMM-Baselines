# CMake generated Testfile for 
# Source directory: /data/lsl/SpGEMM/kokkos-kernels/sparse/unit_test
# Build directory: /data/lsl/SpGEMM/kokkos-kernels/sparse/unit_test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(sparse_serial "/data/lsl/SpGEMM/kokkos-kernels/sparse/unit_test/KokkosKernels_sparse_serial")
set_tests_properties(sparse_serial PROPERTIES  _BACKTRACE_TRIPLES "/data/lsl/SpGEMM/kokkos-kernels/cmake/fake_tribits.cmake;212;add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;194;kokkoskernels_add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;146;kokkoskernels_add_executable_and_test;/data/lsl/SpGEMM/kokkos-kernels/sparse/unit_test/CMakeLists.txt;65;kokkoskernels_add_unit_test;/data/lsl/SpGEMM/kokkos-kernels/sparse/unit_test/CMakeLists.txt;0;")
add_test(blocksparse_serial "/data/lsl/SpGEMM/kokkos-kernels/sparse/unit_test/KokkosKernels_blocksparse_serial")
set_tests_properties(blocksparse_serial PROPERTIES  _BACKTRACE_TRIPLES "/data/lsl/SpGEMM/kokkos-kernels/cmake/fake_tribits.cmake;212;add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;194;kokkoskernels_add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;146;kokkoskernels_add_executable_and_test;/data/lsl/SpGEMM/kokkos-kernels/sparse/unit_test/CMakeLists.txt;69;kokkoskernels_add_unit_test;/data/lsl/SpGEMM/kokkos-kernels/sparse/unit_test/CMakeLists.txt;0;")
add_test(sparse_openmp "/data/lsl/SpGEMM/kokkos-kernels/sparse/unit_test/KokkosKernels_sparse_openmp")
set_tests_properties(sparse_openmp PROPERTIES  _BACKTRACE_TRIPLES "/data/lsl/SpGEMM/kokkos-kernels/cmake/fake_tribits.cmake;212;add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;194;kokkoskernels_add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;146;kokkoskernels_add_executable_and_test;/data/lsl/SpGEMM/kokkos-kernels/sparse/unit_test/CMakeLists.txt;75;kokkoskernels_add_unit_test;/data/lsl/SpGEMM/kokkos-kernels/sparse/unit_test/CMakeLists.txt;0;")
add_test(blocksparse_openmp "/data/lsl/SpGEMM/kokkos-kernels/sparse/unit_test/KokkosKernels_blocksparse_openmp")
set_tests_properties(blocksparse_openmp PROPERTIES  _BACKTRACE_TRIPLES "/data/lsl/SpGEMM/kokkos-kernels/cmake/fake_tribits.cmake;212;add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;194;kokkoskernels_add_test;/data/lsl/SpGEMM/kokkos-kernels/cmake/kokkoskernels_tribits.cmake;146;kokkoskernels_add_executable_and_test;/data/lsl/SpGEMM/kokkos-kernels/sparse/unit_test/CMakeLists.txt;79;kokkoskernels_add_unit_test;/data/lsl/SpGEMM/kokkos-kernels/sparse/unit_test/CMakeLists.txt;0;")
