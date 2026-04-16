#ifndef _TEST_HASPGEMM_
#define _TEST_HASPGEMM_

#include "HASpGEMM.h"
#include "Check_d_csrmultcsr.h"

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
    for(int i = 1;i <= P_num + E_num;i++){
        if(i <= P_num){
            cumulative += gap_P;
            thread_row[i] =  binary_search(sCub,0,bound,cumulative);
        }
        else{
            cumulative += gap_E;
            thread_row[i] = binary_search(sCub,bound,mA,cumulative);
        }
    }
    thread_row[P_num + E_num] = mA;
}

void Test_HA_d_csrmultcsr(char *mat_name,
                       int nA, //n
                       int mA, //m
                       int nB, //m
                       int mB, //n
                       int *csrRowPtrA_1based,
                       int *csrColIdxA_1based,
                       MAT_VAL_TYPE *csrValA,
                       int *csrRowPtrB_1based,
                       int *csrColIdxB_1based,
                       MAT_VAL_TYPE *csrValB,
                       int *nnzCub, double *avg_nnzCub, int *max_nnzCub, int *min_nnzCub,
                       MAT_VAL_TYPE *avg_time, MAT_VAL_TYPE *min_time,
                       MAT_VAL_TYPE *avg_gflops, MAT_VAL_TYPE *max_gflops,
                       MAT_VAL_TYPE *analyze_time,
                       MAT_VAL_TYPE *C_row_variation, int *repeat_value)
{
    int core_num = P_NUM + E_NUM;
    int *thread_row = (int*)malloc((core_num + 1)*sizeof(int));
    memset(thread_row,0,sizeof thread_row);
    struct timeval t1, t2;
    MAT_PTR_TYPE nnzc = 0;
    MAT_PTR_TYPE *csrRowPtrC_1based = (MAT_PTR_TYPE *)malloc((mA + 1) * sizeof(MAT_PTR_TYPE));
    int *csrColIdxC_1based;
    MAT_VAL_TYPE *csrValC;

    MAT_PTR_TYPE nnzA = csrRowPtrA_1based[mA] - 1;
    MAT_PTR_TYPE nnzB = csrRowPtrB_1based[mB] - 1;

    // *repeat_value = 50;
     if(nnzA < 100000)
        *repeat_value = 30;
    else if(nnzA < 1000000)
        *repeat_value = 5;
    else 
        *repeat_value = 1;

    int request;
    double sum_time1 = 0,sum_time2 = 0,sum_time = 0,analyze_min_time = 0;
    *avg_time = 0, *min_time;
    double max_time = 0;
    int *iCub = (int *)malloc(mA * sizeof(int));
    int *sCub = (int *)malloc(mA * sizeof(int));
    int *weight = (int *)malloc(mA * sizeof(int));
    long *sweight = (long *)malloc(mA * sizeof(long));
    int *ca = (int *)malloc(mA * sizeof(int));
    long *sca = (long *)malloc(mA * sizeof(long));
    int max_nnz,min_nnz,avg_nnz;
    

    for(int i = 0;i < *repeat_value;i++){
        memset(iCub, 0, sizeof(int) * mA);
        memset(sCub, 0, sizeof(int) * mA);
        memset(weight, 0, sizeof(int) * mA);
        memset(sweight, 0, sizeof(long) * mA);
    gettimeofday(&t1,NULL);
#pragma omp parallel for
        for(int i = 0;i < mA;i++){
            int loc = -1;
            weight[i] += csrRowPtrA_1based[i + 1] - csrRowPtrA_1based[i];
            int cnt_cache = 1;
            for (int j = csrRowPtrA_1based[i] - csrRowPtrA_1based[0]; j < csrRowPtrA_1based[i + 1] - csrRowPtrA_1based[0]; j++)
            {
                int rowAtop = csrColIdxA_1based[j] - csrRowPtrA_1based[0];
                iCub[i] += csrRowPtrB_1based[rowAtop + 1] - csrRowPtrB_1based[rowAtop];
                for(int k = csrRowPtrB_1based[rowAtop] - csrRowPtrA_1based[0];k < csrRowPtrB_1based[rowAtop + 1] - csrRowPtrA_1based[0];k++){
                    if(loc == -1){
                        loc = k;
                    }
                    if(k / 16 > loc / 16){
                        cnt_cache++;
                        loc = k;
                    }
                }
            }
            weight[i] = iCub[i] + cnt_cache;
        }

        sCub[0] = iCub[0];
        sweight[0] = weight[0];

        for(int i = 1;i < mA;i++){
            sCub[i] = sCub[i - 1] + iCub[i];
            sweight[i] = sweight[i - 1] + weight[i];
        }

    double pro = 0.50;
    if(E_NUM == 16) pro = 0.66;
    HA_parition_weight_row(P_NUM,E_NUM,pro,thread_row,mA,sweight);

    gettimeofday(&t2,NULL);
    analyze_min_time = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
    if(i == 0) *analyze_time = analyze_min_time;
    else *analyze_time = *analyze_time < analyze_min_time ? *analyze_time : analyze_min_time;
    }


    for(int i = 0;i < *repeat_value;i++){
        memset(csrRowPtrC_1based, 0, (mA + 1) * sizeof(MAT_PTR_TYPE));
        request = 1;
        gettimeofday(&t1, NULL);
        haspgemm(&request, &mA, &nA, &nB,
                (MAT_VAL_TYPE *)csrValA, csrColIdxA_1based, csrRowPtrA_1based,
                (MAT_VAL_TYPE *)csrValB, csrColIdxB_1based, csrRowPtrB_1based,
                (MAT_VAL_TYPE *)csrValC, csrColIdxC_1based, csrRowPtrC_1based,
                &nnzc, thread_row,core_num,iCub);
        gettimeofday(&t2, NULL);
        MAT_VAL_TYPE time1_haspgemm = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
        sum_time1 += time1_haspgemm;
        // Make space for arrays
        nnzc = csrRowPtrC_1based[mA] - 1;
        csrColIdxC_1based = (int *)malloc((nnzc) * sizeof(int));
        csrValC = (MAT_VAL_TYPE *)malloc((nnzc) * sizeof(MAT_VAL_TYPE));
        memset(csrColIdxC_1based, 0, (nnzc) * sizeof(int));
        memset(csrValC, 0, (nnzc) * sizeof(MAT_VAL_TYPE));
        request = 2;

        gettimeofday(&t1, NULL);
        haspgemm(&request, &mA, &nA, &nB,
                (MAT_VAL_TYPE *)csrValA, csrColIdxA_1based, csrRowPtrA_1based,
                (MAT_VAL_TYPE *)csrValB, csrColIdxB_1based, csrRowPtrB_1based,
                (MAT_VAL_TYPE *)csrValC, csrColIdxC_1based, csrRowPtrC_1based,
                &nnzc, thread_row,core_num,iCub);
        gettimeofday(&t2, NULL);
        MAT_VAL_TYPE time2_haspgemm = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
        sum_time2 += time2_haspgemm;
        MAT_VAL_TYPE time = time1_haspgemm + time2_haspgemm;
        if (i == 0) // Getting the minimum time
            *min_time = time;
        sum_time += time;
        *min_time = time < *min_time ? time : *min_time;
    }
    free(iCub);
    // Calculate the relevant information of the matrix C
    double *rowC = (double*)malloc(mA * sizeof(double));
    *nnzCub = 0;
    int cub = 0, cubnum = 0;
    *min_nnzCub = (csrRowPtrB_1based[csrColIdxA_1based[0]] - csrRowPtrB_1based[csrColIdxA_1based[0] - 1]) * nnzA;
    *max_nnzCub = 0;
    int m = 0;
    int i = 0;
    int j = 0;
    while (i < nnzA)
    {
        cubnum = 0;
        while (i < csrRowPtrA_1based[m + 1] - 1)
        {
            int rowB = csrColIdxA_1based[i] - 1;
            cub = csrRowPtrB_1based[rowB + 1] - csrRowPtrB_1based[rowB];
            cubnum += cub;
            i++;
        }
        *nnzCub += cubnum;
        rowC[j] = (double)cubnum;
        if (cubnum > *max_nnzCub)
            *max_nnzCub = cubnum;
        if (cubnum < *min_nnzCub)
            *min_nnzCub = cubnum;
        j++;
        m++;
    }
    *avg_nnzCub = (double)*nnzCub / (double)mA;
    *C_row_variation = matrix_variance(rowC, mA);
    free(rowC);
    *avg_time = sum_time / (double)*repeat_value;
    *avg_time += *analyze_time;
    *min_time += *analyze_time;
    *avg_gflops = (2 * (*nnzCub)) / ((*avg_time) * 1000000);
    *max_gflops = ((*nnzCub) * 2) / ((*min_time) * 1000000);

    //Checking correctness
    // int res = check_d_csrmultcsr(1,nA,mA,nB,mB,csrRowPtrA_1based,csrColIdxA_1based,csrValA,csrRowPtrB_1based,csrColIdxB_1based,csrValB,csrRowPtrC_1based,csrColIdxC_1based,csrValC,nnzc);
    // if(res != 0)   printf("Fail!\n");
}
#endif
