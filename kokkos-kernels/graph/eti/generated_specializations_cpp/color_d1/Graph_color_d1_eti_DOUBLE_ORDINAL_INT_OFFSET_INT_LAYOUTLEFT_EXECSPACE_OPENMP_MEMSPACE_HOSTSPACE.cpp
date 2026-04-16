// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project


#define KOKKOSKERNELS_IMPL_COMPILE_LIBRARY true
#include "KokkosKernels_config.h"

#include "KokkosGraph_color_d1_spec.hpp"
namespace KokkosGraph {
namespace Impl {
KOKKOSGRAPH_COLOR_D1_ETI_SPEC_INST(double,int,int,Kokkos::LayoutLeft,Kokkos::OpenMP,Kokkos::HostSpace)
  } //IMPL 
} //Kokkos
