#include <time.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "HASpGEMM.h"
#include "Check_d_csrmultcsr.h"
#include "common.h"
#include "biio.h"
int binary_search(long *csrRowPtr,int left,int right,long nnz_bound){
    int mid = left + right >> 1;
    while(left < right){
        if(csrRowPtr[mid] > nnz_bound) right = mid;
        else left = mid + 1;
        mid = left + right >> 1;
    }
    return mid;
}

void HA_parition_weight_row(int P_num,
                int E_num,
                double propration,
                int *thread_row,
                int mA,
                long *sCub)
{
    long sum_Cub = sCub[mA - 1];
    long sum_P = sum_Cub * propration;
    long sum_E = sum_Cub - sum_P;
    long gap_P = ceil(sum_P / P_num);
    long gap_E = ceil(sum_E / E_num);
    long cumulative = 0;
    int bound = binary_search(sCub,0,mA,sum_P);
    for(int i = 1;i < P_num;i++){
        cumulative += gap_P;
        thread_row[i] =  binary_search(sCub,0,bound,cumulative);
        
    }
    for(int i = P_num;i < P_num + E_num;i++){
        cumulative += gap_E;
        thread_row[i] = binary_search(sCub,bound,mA,cumulative);
    }

    thread_row[P_num + E_num] = mA;
}


int main(int argc, char **argv)
{
    if (argc < 2)
    {
        printf("Run the code by './thmkl_spgemm matrixA.mtx matrixB.mtx'.\n");
        return 0;
    }

    char *filenameA, *filenameB;
    filenameA = argv[1]; // matrixA
    filenameB = argv[1]; // matrixB

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

    // Gets the name of the matrix for recording matrix information
    char *tempA = strtok(filenameA, "/");   
    char *mat_nameA;
    char *mat_group_nameA;
    int t = 0;
    int t_for = 4;
    while (tempA)
    {
        tempA = strtok(NULL, "/");
        if (tempA[strlen(tempA) - 4] == '.')
        {
            mat_nameA = strtok(tempA, ".");
            break;
        }
        if(t == (t_for-2))
        {
            mat_group_nameA = tempA;
        }
        t++;
    }
    char mat_name[100] = "";
    char mat_group_name[100] = "";
    strcat(mat_name, mat_nameA);
    strcat(mat_group_name, mat_group_nameA);

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


    int core_num = 24;
     int *thread_row = (int*)malloc((core_num + 1)*sizeof(int));
    memset(thread_row,0,sizeof thread_row);
    struct timeval t1, t2;
    
    MAT_PTR_TYPE *csrRowPtrC_1based = (MAT_PTR_TYPE *)malloc((mA + 1) * sizeof(MAT_PTR_TYPE));
    int *csrColIdxC_1based;
    MAT_VAL_TYPE *csrValC;

    repeat_value = 1;

    int request;
    double sum_time1 = 0,sum_time2 = 0,sum_time = 0,analyze_min_time = 0;
    avg_time = 0, min_time;
    double max_time = 0;
    int *iCub = (int *)malloc(mA * sizeof(int));
    int *sCub = (int *)malloc(mA * sizeof(int));
    int *weight = (int *)malloc(mA * sizeof(int));
    long *sweight = (long *)malloc(mA * sizeof(long));
    int *ca = (int *)malloc(mA * sizeof(int));
    long *sca = (long *)malloc(mA * sizeof(long));
    int max_nnz,min_nnz,avg_nnz;

    memset(iCub, 0, sizeof(int) * mA);
    memset(sCub, 0, sizeof(int) * mA);
    memset(weight, 0, sizeof(int) * mA);
    memset(sweight, 0, sizeof(long) * mA);
    memset(ca, 0, sizeof(int) * mA);
    memset(sca, 0, sizeof(long) * mA);
    for(int i = 0;i < mA;i++){
             int loc = -1;
            // weight[i] += csrRowPtrA_1based[i + 1] - csrRowPtrA_1based[i];
            int cnt_cache = 1;
            for (int j = csrRowPtrA_1based[i] - csrRowPtrA_1based[0]; j < csrRowPtrA_1based[i + 1] - csrRowPtrA_1based[0]; j++)
            {
                if(loc == -1){
                    loc = j;
                }
                if(j / 16 > loc / 16){
                    cnt_cache++;
                    loc = j;
                }
                int rowAtop = csrColIdxA_1based[j] - csrRowPtrA_1based[0];
                iCub[i] += csrRowPtrB_1based[rowAtop + 1] - csrRowPtrB_1based[rowAtop];
                int lloc = -1;
                int cnt_cache_B = 1;
                for(int k = csrRowPtrB_1based[rowAtop] - csrRowPtrA_1based[0];k < csrRowPtrB_1based[rowAtop + 1] - csrRowPtrA_1based[0];k++){
                    if(lloc == -1){
                        lloc = k;
                    }
                    if(k / 16 > lloc / 16){
                        cnt_cache_B++;
                        lloc = k;
                    }
                }
                cnt_cache += cnt_cache_B;
            }
            weight[i] = iCub[i] + cnt_cache;
            // weight[i] = cnt_cache;
            ca[i] = cnt_cache;
        }
        sCub[0] = iCub[0];
        sweight[0] = weight[0];
        sca[0] = ca[0];

        for(int i = 1;i < mA;i++){
            sCub[i] = sCub[i - 1] + iCub[i];
            sweight[i] = sweight[i - 1] + weight[i];
            sca[i] = sca[i - 1] + ca[i];
            // printf("%d %d\n",i,sCub[i]);
        }

    HA_parition_weight_row(16,8,0.66,thread_row,mA,sweight);

    int gap_ma = mA / 24;
    int cum = 0;
    for(int i = 0; i < 24; i++){
        int cum1 = cum;
        cum += gap_ma;
        if(i == 23) {
            cum = mA - 1;
            printf("%d,%d,%ld,%d,%ld\n",i + 1,sCub[thread_row[i]] / 2,sca[thread_row[i + 1] - 1] - sca[thread_row[i]],sCub[cum] - sCub[cum1],sca[cum] - sca[cum1]);
        }
        else printf("%d,%d,%ld,%d,%ld\n",i + 1,sCub[thread_row[i + 1]] - sCub[thread_row[i]],sca[thread_row[i + 1]] - sca[thread_row[i]],sCub[cum] - sCub[cum1],sca[cum] - sca[cum1]);
    }
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
