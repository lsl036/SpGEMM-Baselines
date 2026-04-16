//------------------------------------------------------------------------------
// GraphBLAS/Demo/Program/spgemm_demo.c: simple SpGEMM (GrB_mxm) benchmark
//------------------------------------------------------------------------------
//
// SuiteSparse:GraphBLAS, Timothy A. Davis, (c) 2017-2025, All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0
//
//------------------------------------------------------------------------------
//
// Usage:
//
//      spgemm_demo Afile Bfile [ntrials]
//
// Afile, Bfile: matrices stored either as:
//      (1) Matrix Market format (.mtx), with 1-based indices, or
//      (2) plain text triplets "i j x" per line, with 0-based indices.
// For Matrix Market "pattern" matrices, all values are taken as 1.0.
//
// The program reads sparse matrices A and B in double precision, checks that
// their inner dimensions are compatible, and then computes C = A*B using
// GrB_mxm on the standard plus-times semiring.  The multiplication is
// repeated ntrials times (default: 5) and the best time is reported.
//
// This demo is meant only as a simple SpGEMM performance test.
//

#include "graphblas_demos.h"

// FREE_ALL is used by the OK(...) macro from graphblas_demos.h.
#undef FREE_ALL
#define FREE_ALL                       \
    GrB_Matrix_free (&A) ;             \
    GrB_Matrix_free (&B) ;             \
    GrB_Matrix_free (&C) ;

//------------------------------------------------------------------------------
// read a Matrix Market file (real or pattern) into an FP64 matrix
//------------------------------------------------------------------------------

static GrB_Info read_matrix_market_fp64   // read Matrix Market into GrB_FP64
(
    GrB_Matrix *A_output,               // matrix to create
    FILE *f,                            // already-open file
    bool pr                             // if true, print basic info
)
{
    if (A_output == NULL || f == NULL)
    {
        return (GrB_NULL_POINTER) ;
    }

    char line [256] ;

    // read header line
    if (fgets (line, sizeof (line), f) == NULL)
    {
        return (GrB_INVALID_VALUE) ;
    }

    // expect "%%MatrixMarket matrix coordinate ..."
    if (strncmp (line, "%%MatrixMarket", 14) != 0)
    {
        return (GrB_INVALID_VALUE) ;
    }

    bool is_pattern = (strstr (line, "pattern") != NULL) ;

    // skip comment lines starting with '%'
    long pos = ftell (f) ;
    int c = fgetc (f) ;
    while (c == '%')
    {
        if (fgets (line, sizeof (line), f) == NULL)
        {
            return (GrB_INVALID_VALUE) ;
        }
        pos = ftell (f) ;
        c = fgetc (f) ;
    }
    // rewind one character to start of dimensions line
    if (c != EOF)
    {
        fseek (f, pos, SEEK_SET) ;
    }

    // read dimensions: nrows ncols nnz
    int64_t nrows, ncols, ntuples ;
    if (fscanf (f, "%ld %ld %ld", &nrows, &ncols, &ntuples) != 3)
    {
        return (GrB_INVALID_VALUE) ;
    }

    if (pr)
    {
        printf ("MatrixMarket: nrows %.16g ncols %.16g ntuples %.16g%s\n",
            (double) nrows, (double) ncols, (double) ntuples,
            is_pattern ? " (pattern)" : "") ;
    }

    if (ntuples < 0 || nrows <= 0 || ncols <= 0)
    {
        return (GrB_INVALID_VALUE) ;
    }

    GrB_Index nnz = (GrB_Index) ntuples ;

    // allocate arrays
    GrB_Index *I = (GrB_Index *) malloc (nnz * sizeof (GrB_Index)) ;
    GrB_Index *J = (GrB_Index *) malloc (nnz * sizeof (GrB_Index)) ;
    double    *X = (double    *) malloc (nnz * sizeof (double   )) ;

    if (I == NULL || J == NULL || X == NULL)
    {
        if (I != NULL) free (I) ;
        if (J != NULL) free (J) ;
        if (X != NULL) free (X) ;
        return (GrB_OUT_OF_MEMORY) ;
    }

    // read entries, 1-based indices in the file
    for (GrB_Index k = 0 ; k < nnz ; k++)
    {
        int64_t i2, j2 ;
        double x = 1.0 ;

        if (is_pattern)
        {
            if (fscanf (f, "%ld %ld", &i2, &j2) != 2)
            {
                free (I) ; free (J) ; free (X) ;
                return (GrB_INVALID_VALUE) ;
            }
        }
        else
        {
            if (fscanf (f, "%ld %ld %lg", &i2, &j2, &x) != 3)
            {
                free (I) ; free (J) ; free (X) ;
                return (GrB_INVALID_VALUE) ;
            }
        }

        // convert to 0-based for GraphBLAS
        i2-- ;
        j2-- ;
        if (i2 < 0 || j2 < 0 || i2 >= nrows || j2 >= ncols)
        {
            free (I) ; free (J) ; free (X) ;
            return (GrB_INVALID_INDEX) ;
        }

        I [k] = (GrB_Index) i2 ;
        J [k] = (GrB_Index) j2 ;
        X [k] = x ;
    }

    // build matrix
    GrB_Matrix A = NULL ;
    GrB_Info info = GrB_Matrix_new (&A, GrB_FP64,
        (GrB_Index) nrows, (GrB_Index) ncols) ;
    if (info != GrB_SUCCESS)
    {
        free (I) ; free (J) ; free (X) ;
        return (info) ;
    }

    info = GrB_Matrix_build_FP64 (A, I, J, X, nnz, GrB_PLUS_FP64) ;
    free (I) ; free (J) ; free (X) ;
    if (info != GrB_SUCCESS)
    {
        GrB_Matrix_free (&A) ;
        return (info) ;
    }

    *A_output = A ;
    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// read plain triplets "i j x" (0-based) into an FP64 matrix
//------------------------------------------------------------------------------

static GrB_Info read_triplet_fp64          // read triplet file into GrB_FP64
(
    GrB_Matrix *A_output,                 // matrix to create
    FILE *f,                              // already-open file
    bool one_based,                       // if true, convert 1-based to 0-based
    bool pr                               // if true, print basic info
)
{
    if (A_output == NULL || f == NULL)
    {
        return (GrB_NULL_POINTER) ;
    }

    int64_t len = 256 ;
    int64_t ntuples = 0 ;
    double x ;

    GrB_Index *I = (GrB_Index *) malloc (len * sizeof (GrB_Index)) ;
    GrB_Index *J = (GrB_Index *) malloc (len * sizeof (GrB_Index)) ;
    double    *X = (double    *) malloc (len * sizeof (double   )) ;

    if (I == NULL || J == NULL || X == NULL)
    {
        if (I != NULL) free (I) ;
        if (J != NULL) free (J) ;
        if (X != NULL) free (X) ;
        return (GrB_OUT_OF_MEMORY) ;
    }

    // read triplets until EOF
    double i2, j2 ;
    while (fscanf (f, "%lg %lg %lg", &i2, &j2, &x) == 3)
    {
        int64_t i = (int64_t) i2 ;
        int64_t j = (int64_t) j2 ;

        if (one_based)
        {
            i-- ;
            j-- ;
        }

        if (i < 0 || j < 0)
        {
            free (I) ; free (J) ; free (X) ;
            return (GrB_INVALID_INDEX) ;
        }

        if (ntuples >= len)
        {
            int64_t newlen = len * 2 ;
            GrB_Index *I2 = (GrB_Index *) realloc (I,
                newlen * sizeof (GrB_Index)) ;
            GrB_Index *J2 = (GrB_Index *) realloc (J,
                newlen * sizeof (GrB_Index)) ;
            double    *X2 = (double    *) realloc (X,
                newlen * sizeof (double   )) ;
            if (I2 == NULL || J2 == NULL || X2 == NULL)
            {
                if (I2 != NULL) I = I2 ;
                if (J2 != NULL) J = J2 ;
                if (X2 != NULL) X = X2 ;
                free (I) ; free (J) ; free (X) ;
                return (GrB_OUT_OF_MEMORY) ;
            }
            I = I2 ;
            J = J2 ;
            X = X2 ;
            len = newlen ;
        }

        I [ntuples] = (GrB_Index) i ;
        J [ntuples] = (GrB_Index) j ;
        X [ntuples] = x ;
        ntuples++ ;
    }

    if (ntuples == 0)
    {
        free (I) ; free (J) ; free (X) ;
        return (GrB_INVALID_VALUE) ;
    }

    if (pr) printf ("ntuples: %.16g\n", (double) ntuples) ;

    // determine matrix dimensions
    int64_t nrows = 0 ;
    int64_t ncols = 0 ;
    for (int64_t k = 0 ; k < ntuples ; k++)
    {
        if ((int64_t) I [k] > nrows) nrows = (int64_t) I [k] ;
        if ((int64_t) J [k] > ncols) ncols = (int64_t) J [k] ;
    }
    nrows++ ;
    ncols++ ;

    if (pr) printf ("nrows %.16g ncols %.16g\n",
        (double) nrows, (double) ncols) ;

    GrB_Matrix A = NULL ;
    GrB_Info info = GrB_Matrix_new (&A, GrB_FP64,
        (GrB_Index) nrows, (GrB_Index) ncols) ;
    if (info != GrB_SUCCESS)
    {
        free (I) ; free (J) ; free (X) ;
        return (info) ;
    }

    info = GrB_Matrix_build_FP64 (A, I, J, X, (GrB_Index) ntuples,
        GrB_PLUS_FP64) ;
    free (I) ; free (J) ; free (X) ;
    if (info != GrB_SUCCESS)
    {
        GrB_Matrix_free (&A) ;
        return (info) ;
    }

    *A_output = A ;
    return (GrB_SUCCESS) ;
}

//------------------------------------------------------------------------------
// auto-detect format and read: Matrix Market or plain triplet
//------------------------------------------------------------------------------

static GrB_Info read_matrix_auto
(
    GrB_Matrix *A_output,
    FILE *f,
    bool pr
)
{
    if (A_output == NULL || f == NULL)
    {
        return (GrB_NULL_POINTER) ;
    }

    // peek first non-whitespace character
    int c = fgetc (f) ;
    while (c == ' ' || c == '\t' || c == '\n' || c == '\r')
    {
        c = fgetc (f) ;
    }
    if (c == EOF)
    {
        return (GrB_INVALID_VALUE) ;
    }
    ungetc (c, f) ;

    if (c == '%')
    {
        // likely Matrix Market
        return (read_matrix_market_fp64 (A_output, f, pr)) ;
    }
    else
    {
        // plain triplets "i j x" (assume 0-based)
        return (read_triplet_fp64 (A_output, f, false, pr)) ;
    }
}

int main (int argc, char **argv)
{
    GrB_Info info ;
    GrB_Matrix A = NULL, B = NULL, C = NULL ;

    //--------------------------------------------------------------------------
    // check command line
    //--------------------------------------------------------------------------

    if (argc < 3)
    {
        fprintf (stderr,
            "usage: %s Afile Bfile [ntrials]\n"
            "  Afile, Bfile: matrices in triplet form (i j x), 0-based indices\n"
            "  ntrials: number of repeated SpGEMM runs (default: 5)\n",
            argv [0]) ;
        return (1) ;
    }

    int ntrials = 5 ;
    if (argc > 3)
    {
        ntrials = atoi (argv [3]) ;
        if (ntrials <= 0) ntrials = 5 ;
    }

    //--------------------------------------------------------------------------
    // start GraphBLAS
    //--------------------------------------------------------------------------

    OK (GrB_init (GrB_NONBLOCKING)) ;

    int32_t nthreads ;
    OK (GrB_Global_get_INT32 (GrB_GLOBAL, &nthreads, GxB_NTHREADS)) ;
    fprintf (stderr, "spgemm_demo: nthreads %d\n", nthreads) ;

    //--------------------------------------------------------------------------
    // read matrices A and B from files (Matrix Market or triplet)
    //--------------------------------------------------------------------------

    FILE *fa = fopen (argv [1], "r") ;
    if (fa == NULL)
    {
        fprintf (stderr, "spgemm_demo: cannot open file %s\n", argv [1]) ;
        GrB_finalize ( ) ;
        return (1) ;
    }

    FILE *fb = fopen (argv [2], "r") ;
    if (fb == NULL)
    {
        fprintf (stderr, "spgemm_demo: cannot open file %s\n", argv [2]) ;
        fclose (fa) ;
        GrB_finalize ( ) ;
        return (1) ;
    }

    // read A and B, auto-detecting the file format
    OK (read_matrix_auto (&A, fa, true)) ;
    OK (read_matrix_auto (&B, fb, true)) ;

    fclose (fa) ;
    fclose (fb) ;

    GrB_Index nrowsA, ncolsA, nrowsB, ncolsB ;
    OK (GrB_Matrix_nrows (&nrowsA, A)) ;
    OK (GrB_Matrix_ncols (&ncolsA, A)) ;
    OK (GrB_Matrix_nrows (&nrowsB, B)) ;
    OK (GrB_Matrix_ncols (&ncolsB, B)) ;

    if (ncolsA != nrowsB)
    {
        fprintf (stderr,
            "spgemm_demo: dimension mismatch, A is %" PRIu64 " x %" PRIu64
            ", B is %" PRIu64 " x %" PRIu64 "\n",
            (uint64_t) nrowsA, (uint64_t) ncolsA,
            (uint64_t) nrowsB, (uint64_t) ncolsB) ;
        FREE_ALL ;
        GrB_finalize ( ) ;
        return (1) ;
    }

    //--------------------------------------------------------------------------
    // prepare for SpGEMM
    //--------------------------------------------------------------------------

    // A and B are FP64, use standard plus-times semiring
    GrB_Semiring semiring = GrB_PLUS_TIMES_SEMIRING_FP64 ;

    double best_time = 1e300 ;
    double total_time = 0 ;

    printf ("A: %" PRIu64 " x %" PRIu64 ", B: %" PRIu64 " x %" PRIu64 "\n",
        (uint64_t) nrowsA, (uint64_t) ncolsA,
        (uint64_t) nrowsB, (uint64_t) ncolsB) ;
    printf ("SpGEMM trials: %d\n", ntrials) ;

    //--------------------------------------------------------------------------
    // run SpGEMM multiple times
    //--------------------------------------------------------------------------

    for (int trial = 0 ; trial < ntrials ; trial++)
    {
        GrB_Matrix_free (&C) ;
        OK (GrB_Matrix_new (&C, GrB_FP64, nrowsA, ncolsB)) ;

        double t1 = WALLCLOCK ;

        OK (GrB_mxm (C, NULL, NULL, semiring, A, B, NULL)) ;

        // ensure computation is finished for accurate timing
        OK (GrB_Matrix_wait (C, GrB_MATERIALIZE)) ;

        double t2 = WALLCLOCK ;
        double dt = t2 - t1 ;
        total_time += dt ;
        if (dt < best_time) best_time = dt ;

        GrB_Index nnzC ;
        OK (GrB_Matrix_nvals (&nnzC, C)) ;

        printf ("trial %d: time %g sec, nnz(C) = %" PRIu64 "\n",
            trial, dt, (uint64_t) nnzC) ;
    }

    printf ("best time: %g sec\n", best_time) ;
    printf ("avg  time: %g sec\n", total_time / ntrials) ;

    //--------------------------------------------------------------------------
    // clean up and finish
    //--------------------------------------------------------------------------

    FREE_ALL ;
    GrB_finalize ( ) ;
    return (0) ;
}

