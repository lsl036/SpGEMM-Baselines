// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#define KOKKOSKERNELS_IMPL_COMPILE_LIBRARY true
#include "KokkosBatched_HostLevel_Gemm.hpp"
namespace KokkosBatched {
namespace Impl {
using KokkosBlas::Trans;
KOKKOSBATCHED_GEMM_T_NT_BLR_ETI_SPEC_INST(double,Kokkos::LayoutLeft,Kokkos::OpenMP,Kokkos::HostSpace)
}  // namespace Impl
}  // namespace KokkosBatched
