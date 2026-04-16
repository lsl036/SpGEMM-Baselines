// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project


#define KOKKOSKERNELS_IMPL_COMPILE_LIBRARY true
#include "KokkosKernels_config.h"
#include "KokkosBlas1_abs_spec.hpp"

namespace KokkosBlas {
namespace Impl {

KOKKOSBLAS1_ABS_ETI_SPEC_INST(double,Kokkos::LayoutLeft,Kokkos::OpenMP,Kokkos::HostSpace)

} // Impl
} // KokkosBlas

