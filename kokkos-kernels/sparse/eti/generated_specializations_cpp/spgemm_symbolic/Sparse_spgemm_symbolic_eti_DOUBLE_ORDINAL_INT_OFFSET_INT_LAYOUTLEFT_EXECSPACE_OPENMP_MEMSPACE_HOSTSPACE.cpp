// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project


#define KOKKOSKERNELS_IMPL_COMPILE_LIBRARY true
#include "KokkosKernels_config.h"

#include "KokkosSparse_spgemm_symbolic_spec.hpp"
namespace KokkosSparse {
namespace Impl {
KOKKOSSPARSE_SPGEMM_SYMBOLIC_ETI_SPEC_INST(double,int,int,Kokkos::LayoutLeft,Kokkos::OpenMP,Kokkos::HostSpace)
  } //IMPL 
} //Kokkos
