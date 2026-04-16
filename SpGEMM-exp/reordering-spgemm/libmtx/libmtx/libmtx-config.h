/* libmtx/libmtx-config.h.  Generated from libmtx-config.h.in by configure.  */
/* This file is part of Libmtx.
 *
 * Copyright (C) 2023 James D. Trotter
 *
 * Libmtx is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * Libmtx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Libmtx.  If not, see
 * <https://www.gnu.org/licenses/>.
 *
 * Authors: James D. Trotter <james@simula.no>
 * Last modified: 2023-06-09
 *
 * Build configuration
 */

#ifndef LIBMTX_LIBMTX_CONFIG_H
#define LIBMTX_LIBMTX_CONFIG_H

/* Define if you have a BLAS library. */
#define LIBMTX_HAVE_BLAS 1

/* Define if you have an MPI library. */
#define LIBMTX_HAVE_MPI 1

/* Define if you have `z' library (-lz) */
#define LIBMTX_HAVE_LIBZ 1

/* Define if you have libpng */
#define LIBMTX_HAVE_LIBPNG 1

/* Define if you have METIS */
#define LIBMTX_HAVE_METIS 1

/* Define if you have SCOTCH */
/* #undef LIBMTX_HAVE_SCOTCH */

/* Define if you have OpenBLAS. */
#define LIBMTX_HAVE_OPENBLAS 1

/* Define if you have OpenMP. */
#define LIBMTX_HAVE_OPENMP 1

/* Define to 1 or 0, depending whether the compiler supports simple visibility
   declarations. */
#define HAVE_VISIBILITY 1

/* symbol visibility */
#if HAVE_VISIBILITY && LIBMTX_API_EXPORT
#define LIBMTX_API __attribute__((__visibility__("default")))
#else
#define LIBMTX_API
#endif

#endif
