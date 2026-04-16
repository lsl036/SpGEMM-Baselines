// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project


#define KOKKOSKERNELS_IMPL_COMPILE_LIBRARY true
#include "KokkosKernels_config.h"
#include "KokkosLapack_trtri_spec.hpp"

namespace KokkosLapack {
namespace Impl {
KOKKOSLAPACK_TRTRI_ETI_SPEC_INST(double,Kokkos::LayoutLeft,Kokkos::OpenMP,Kokkos::HostSpace)
  } //IMPL 
} //Kokkos
