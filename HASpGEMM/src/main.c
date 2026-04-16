#include <time.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include "common.h"
#include "biio.h"
#include "Test_HASpGEMM.h"

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        printf("Run the code by './haspgemm12 matrixA.mtx matrixB.mtx'.\n");
        return 0;
    }

    char *filenameA, *filenameB;
    filenameA = argv[1]; // matrixA
    filenameB = (argc >= 3) ? argv[2] : argv[1]; // matrixB (default: same as A)

    int mA;
    int nA;
    int nnzA;
    int isSymmetricA;
    MAT_PTR_TYPE *csrRowPtrA_0based;
    int *csrColIdxA_0based;
    MAT_VAL_TYPE *csrValA;
    // get matrix A (support cbd file or mtx file)
    read_Dmatrix_32(&mA, &nA, &nnzA, &csrRowPtrA_0based, &csrColIdxA_0based, &csrValA, &isSymmetricA, filenameA);
    srand(1);


    int mB;
    int nB;
    int nnzB;
    int isSymmetricB;
    MAT_PTR_TYPE *csrRowPtrB_0based;
    int *csrColIdxB_0based;
    MAT_VAL_TYPE *csrValB;

    mB = nA;
    nB = mA;
    nnzB = nnzA;
    csrRowPtrB_0based = (MAT_PTR_TYPE *)malloc(sizeof(MAT_PTR_TYPE) * (mB + 1));
    csrColIdxB_0based = (int *)malloc(sizeof(int) * nnzB);
    csrValB = (MAT_VAL_TYPE *)malloc(sizeof(MAT_VAL_TYPE) * nnzB);
    // Get matrix B 
    matrix_transposition(mA, nA, nnzA,
                         csrRowPtrA_0based, csrColIdxA_0based, csrValA,
                         csrColIdxB_0based, csrRowPtrB_0based, csrValB);

    float *csrValBs = (float *)malloc(sizeof(float) * nnzB);
    for (int i = 0; i < nnzB; i++)
        csrValBs[i] = csrValB[i];

    // Gets the name of the matrix for recording matrix information (use copy so argv[1] is not modified)
    char path_buf[1024];
    size_t path_len = strlen(filenameA);
    if (path_len >= sizeof(path_buf)) path_len = sizeof(path_buf) - 1;
    memcpy(path_buf, filenameA, path_len);
    path_buf[path_len] = '\0';

    char *mat_nameA = "matrix";
    char *mat_group_nameA = "";
    char *last_part = NULL;
    char *prev_part = NULL;
    char *tok = strtok(path_buf, "/");
    while (tok != NULL)
    {
        prev_part = last_part;
        last_part = tok;
        tok = strtok(NULL, "/");
    }
    if (last_part != NULL)
    {
        mat_nameA = last_part;
        mat_group_nameA = (prev_part != NULL) ? prev_part : last_part;
    }
    char mat_name[100] = "";
    char mat_group_name[100] = "";
    if (mat_nameA != NULL)
    {
        size_t n = 0;
        while (n < sizeof(mat_name) - 1 && mat_nameA[n] != '\0' && mat_nameA[n] != '.') { mat_name[n] = mat_nameA[n]; n++; }
        mat_name[n] = '\0';
    }
    if (mat_group_nameA != NULL)
        strncat(mat_group_name, mat_group_nameA, sizeof(mat_group_name) - 1);

    // create matrixA-1based csr
    MAT_PTR_TYPE *csrRowPtrA_1based = (MAT_PTR_TYPE *)malloc((mA + 1) * sizeof(MAT_PTR_TYPE));
    int *csrColIdxA_1based = (int *)malloc(nnzA * sizeof(int));
#pragma omp parallel for
    for (MAT_PTR_TYPE i = 0; i < nnzA; i++)
        csrColIdxA_1based[i] = csrColIdxA_0based[i] + 1;
#pragma omp parallel for
    for (MAT_PTR_TYPE i = 0; i < mA + 1; i++)
        csrRowPtrA_1based[i] = csrRowPtrA_0based[i] + 1;

    // create matrixB-1based csr
    MAT_PTR_TYPE *csrRowPtrB_1based = (MAT_PTR_TYPE *)malloc((mB + 1) * sizeof(MAT_PTR_TYPE));
    int *csrColIdxB_1based = (int *)malloc(nnzB * sizeof(int));
#pragma omp parallel for
    for (MAT_PTR_TYPE i = 0; i < nnzA; i++)
        csrColIdxB_1based[i] = csrColIdxB_0based[i] + 1;
#pragma omp parallel for
    for (MAT_PTR_TYPE i = 0; i < mB + 1; i++)
        csrRowPtrB_1based[i] = csrRowPtrB_0based[i] + 1;


    int nnzCub;
    double avg_nnzCub;
    int max_nnzCub, min_nnzCub;
    double avg_time, min_time, avg_gflops, max_gflops, analyze_time;
    int repeat_value;
    double C_row_variation;

    // test 
    Test_HA_d_csrmultcsr(mat_name, nA, mA, nB, mB, csrRowPtrA_1based, csrColIdxA_1based, csrValA, csrRowPtrB_1based, csrColIdxB_1based, csrValB, &nnzCub, &avg_nnzCub, &max_nnzCub, &min_nnzCub, &avg_time, &min_time, &avg_gflops, &max_gflops, &analyze_time, &C_row_variation, &repeat_value);

    printf("HASpGEMM,%s,%d,%d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf\n",mat_name,nnzA,mA,avg_nnzCub,avg_gflops,max_gflops,avg_time,analyze_time);

    free(csrRowPtrA_0based);
    free(csrColIdxA_0based);
    free(csrValA);

    free(csrRowPtrB_0based);
    free(csrColIdxB_0based);
    free(csrValB);
    free(csrValBs);

    free(csrRowPtrA_1based);
    free(csrColIdxA_1based);
    free(csrRowPtrB_1based);
    free(csrColIdxB_1based);

    return 0;
}
