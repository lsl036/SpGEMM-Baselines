// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#define KOKKOSKERNELS_IMPL_COMPILE_LIBRARY true
#include "KokkosKernels_config.h"
#include "KokkosBlas2_syr2_spec.hpp"

namespace KokkosBlas {
namespace Impl {
KOKKOSBLAS2_SYR2_ETI_SPEC_INST(double,Kokkos::LayoutLeft,Kokkos::OpenMP,Kokkos::HostSpace)
} //IMPL
} //Kokkos
