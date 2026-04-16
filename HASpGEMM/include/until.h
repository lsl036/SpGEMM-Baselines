#ifndef UNTIL_H_
#define UNTIL_H_
// #include "../src/mmio.h"
// #include "../src/mmio_highlevel.h"
#include "common.h"
#include <omp.h>
#include <stdio.h>
#include <math.h>
// #ifdef ARM_CPU
// #include <arm_neon.h>
// #else
// #include <immintrin.h>
// #include <emmintrin.h>
// #endif
#define cur_type double

#define MKL_INT int
#define CBLAS_INDEX size_t

void swap_key(int *a, int *b)
{
    int tmp = *a;
    if (a != NULL && b != NULL)
    {
        *a = *b;
    }
    if (b != NULL)
    {
        *b = tmp;
    }
}
void swap_val(MAT_VAL_TYPE *a, MAT_VAL_TYPE *b)
{
    MAT_VAL_TYPE tmp = *a;
    if (a != NULL && b != NULL)
    {
        *a = *b;
    }
    if (b != NULL)
    {
        *b = tmp;
    }
}
int partition_key(int *key, int length, int pivot_index, int mid)
{
    int i = 0;
    int small_length = pivot_index;

    int pivot = key[mid];
    if (key != NULL)
    {
        swap_key(&key[pivot_index], &key[pivot_index + (length - 1)]);
    }

    for (; i < length; i++)
    {
        if (key != NULL)
            if (key[pivot_index + i] < pivot)
            {
                swap_key(&key[pivot_index + i], &key[small_length]);
                small_length++;
            }
    }

    if (key != NULL)
    {
        swap_key(&key[pivot_index + length - 1], &key[small_length]);
    }

    return small_length;
}

void insert_sort_key(int *key, int length)
{
    for (int i = 1; i < length; i++)
    {
        int tmp_key = key[i];
        int j = i - 1;
        while ((j >= 0) && (key[j] > tmp_key))
        {
            key[j + 1] = key[j];
            j--;
        }
        key[j + 1] = tmp_key;
    }
}
void quick_sort_key1(int *key, int length)
{
    if (length == 0 || length == 1)
        return;

    int sorted = 1;

    // get mid of three
    int first = 0;
    int last = length - 1;
    int mid = first + ((last - last) >> 1);
    // int pivot = mid;
    if (key[mid] > key[first])
    {
        swap_key(&key[mid], &key[first]);
    }
    if (key[first] > key[last])
    {
        swap_key(&key[first], &key[last]);
    }
    if (key[mid] > key[last])
    {
        swap_key(&key[mid], &key[last]);
    }

    for (int i = 1; i < length; i++)
    {
        if (key[i] < key[i - 1])
        {
            sorted = 0;
            break;
        }
    }

    if (!sorted)
    {
        if (length > 64)
        {
            int small_length = partition_key(key, length, 0, mid);
            quick_sort_key1(key, small_length);
            quick_sort_key1(&key[small_length + 1], length - small_length - 1);
        }
        else
        {
            insert_sort_key(key, length);
        }
    }
}
void scan_par(MAT_PTR_TYPE *a, int SIZE)
{
    int num_threads = 0;
    int seq_num = 0;
    int i, j, tid;
#pragma omp parallel private(tid, j) shared(a)
    {
        tid = omp_get_thread_num();
        num_threads = omp_get_num_threads();
        seq_num = SIZE / num_threads;

        for (j = 0; j < seq_num - 1; j++)
        {
            a[tid * seq_num + j + 1] += a[tid * seq_num + j];
        }
    }

    for (i = 2; i < num_threads + 1; i++)
    {
        a[i * seq_num - 1] += a[(i - 1) * seq_num - 1];
    }

#pragma omp parallel private(tid, j) shared(a)
    {
        tid = omp_get_thread_num();
        if (tid != 0)
        {
            for (j = 0; j < seq_num - 1; j++)
            {
                a[tid * seq_num + j] += a[tid * seq_num - 1];
            }
        }
    }
    if (SIZE % num_threads != 0)
    {
        for (i = num_threads * seq_num; i < SIZE; i++)
        {
            a[i] += a[i - 1];
        }
    }
}
void exclusive_scan(MAT_PTR_TYPE *input, int length)
{
    if (length == 0 || length == 1)
        return;

    MAT_PTR_TYPE old_val, new_val;

    old_val = input[0];
    input[0] = 0;
    for (int i = 1; i < length; i++)
    {
        new_val = input[i];
        input[i] = old_val + input[i - 1];
        old_val = new_val;
    }
}
int partition_key_val_pair(int *key, MAT_VAL_TYPE *val, int length, int pivot_index, int mid)
{
    int i = 0;
    int small_length = pivot_index;

    int pivot = key[pivot_index];
    if (val != NULL && key != NULL)
    {
        swap_key(&key[pivot_index], &key[pivot_index + (length - 1)]);
        swap_val(&val[pivot_index], &val[pivot_index + (length - 1)]);
    }

    for (; i < length; i++)
    {
        if (key != NULL)
        {
            if (key[pivot_index + i] < pivot)
            {
                swap_key(&key[pivot_index + i], &key[small_length]);
                if (val != NULL)
                {
                    swap_val(&val[pivot_index + i], &val[small_length]);
                }
                small_length++;
            }
        }
    }

    if (key != NULL)
    {
        swap_key(&key[pivot_index + length - 1], &key[small_length]);
    }
    if (val != NULL)
    {
        swap_val(&val[pivot_index + length - 1], &val[small_length]);
    }

    return small_length;
}
void insert_sort_key_val_pair(int *key, MAT_VAL_TYPE *val, int length)
{
    for (int i = 1; i < length; i++)
    {
        int tmp_key = key[i];
        MAT_VAL_TYPE tmp_val = val[i];
        int j = i - 1;
        while ((j >= 0) && (key[j] > tmp_key))
        {
            key[j + 1] = key[j];
            val[j + 1] = val[j];
            j--;
        }
        key[j + 1] = tmp_key;
        val[j + 1] = tmp_val;
    }
}

void merge_sort_key_val_pair1(int *key, MAT_VAL_TYPE *val, int length)
{
    if (length <= 1)
        return;

    int mid = length / 2;
    merge_sort_key_val_pair1(key, val, mid);
    merge_sort_key_val_pair1(key + mid, val + mid, length - mid);

    int *tmp_key = (int*)malloc(length * sizeof(int));
    MAT_VAL_TYPE *tmp_val = (MAT_VAL_TYPE*)malloc(length * sizeof(MAT_VAL_TYPE));

    int i = 0, j = mid, k = 0;
    while (i < mid && j < length)
    {
        if (key[i] <= key[j])
        {
            tmp_key[k] = key[i];
            tmp_val[k] = val[i];
            i++;
        }
        else
        {
            tmp_key[k] = key[j];
            tmp_val[k] = val[j];
            j++;
        }
        k++;
    }

    while (i < mid)
    {
        tmp_key[k] = key[i];
        tmp_val[k] = val[i];
        i++;
        k++;
    }

    while (j < length)
    {
        tmp_key[k] = key[j];
        tmp_val[k] = val[j];
        j++;
        k++;
    }

    for (int i = 0; i < length; i++)
    {
        key[i] = tmp_key[i];
        val[i] = tmp_val[i];
    }

    free(tmp_key);
    free(tmp_val);
}

void quick_sort_key_val_pair1(int *key, MAT_VAL_TYPE *val, int length)
{
    if (length == 0 || length == 1)
        return;

    int sorted = 1;

    // get mid of three
    int first = 0;
    int last = length - 1;
    int mid = first + ((last - last) >> 1);
    // int pivot = mid;
    if (key[mid] > key[first])
    {
        swap_key(&key[mid], &key[first]);
        swap_val(&val[mid], &val[first]);
    }
    if (key[first] > key[last])
    {
        swap_key(&key[first], &key[last]);
        swap_val(&val[first], &val[last]);
    }
    if (key[mid] > key[last])
    {
        swap_key(&key[mid], &key[last]);
        swap_val(&val[mid], &val[last]);
    }

    for (int i = 1; i < length; i++)
    {
        if (key[i] < key[i - 1])
        {
            sorted = 0;
            break;
        }
    }

    if (!sorted)
    {
        if (length > 64)
        {
            int small_length = partition_key_val_pair(key, val, length, 0, mid);
            quick_sort_key_val_pair1(key, val, small_length);
            quick_sort_key_val_pair1(&key[small_length + 1], &val[small_length + 1], length - small_length - 1);
        }
        else
        {
            insert_sort_key_val_pair(key, val, length);
        }
    }
}
void segmented_sum(MAT_VAL_TYPE *input, int *bit_flag, int length)
{
    if (length == 0 || length == 1)
        return;

    for (int i = 0; i < length; i++)
    {
        if (bit_flag[i])
        {
            int j = i + 1;
            while (!bit_flag[j] && j < length)
            {
                input[i] += input[j];
                j++;
            }
        }
    }

}

void merge_sort_key1(int *key, int length) {
    if (length <= 1) {
        return;
    }

    int mid = length / 2;

    // 递归排序左半部分
    merge_sort_key1(key, mid);

    // 递归排序右半部分
    merge_sort_key1(key + mid, length - mid);

    // 合并左右两个有序数组
    int *temp = (int*)malloc(length * sizeof(int));
    int i = 0, j = mid, k = 0;
    while (i < mid && j < length) {
        if (key[i] <= key[j]) {
            temp[k++] = key[i++];
        } else {
            temp[k++] = key[j++];
        }
    }
    while (i < mid) {
        temp[k++] = key[i++];
    }
    while (j < length) {
        temp[k++] = key[j++];
    }
    memcpy(key, temp, length * sizeof(int));
    free(temp);
}

void swap_val_s(float *a, float *b)
{
    float tmp = *a;
    if (a != NULL && b != NULL)
    {
        *a = *b;
    }
    if (b != NULL)
    {
        *b = tmp;
    }
}
void insert_sort_key_val_pair_s(int *key, float *val, int length)
{
    for (int i = 1; i < length; i++)
    {
        int tmp_key = key[i];
        float tmp_val = val[i];
        int j = i - 1;
        while ((j >= 0) && (key[j] > tmp_key))
        {
            key[j + 1] = key[j];
            val[j + 1] = val[j];
            j--;
        }
        key[j + 1] = tmp_key;
        val[j + 1] = tmp_val;
    }
}
void matrix_transposition(const int m,
                          const int n,
                          const MAT_PTR_TYPE nnz,
                          const MAT_PTR_TYPE *csrRowPtr,
                          const int *csrColIdx,
                          const MAT_VAL_TYPE *csrVal,
                          int *cscRowIdx,
                          MAT_PTR_TYPE *cscColPtr,
                          MAT_VAL_TYPE *cscVal)
{
    // histogram in column pointer
    memset(cscColPtr, 0, sizeof(MAT_PTR_TYPE) * (n + 1));
    for (MAT_PTR_TYPE i = 0; i < nnz; i++)
    {
        if (cscColPtr != NULL && csrColIdx != NULL && csrRowPtr != NULL)
        {
            cscColPtr[csrColIdx[i] - csrRowPtr[0]]++;
        }
    }

    // prefix-sum scan to get the column pointer
    exclusive_scan(cscColPtr, n + 1);

    MAT_PTR_TYPE *cscColIncr = (MAT_PTR_TYPE *)malloc(sizeof(MAT_PTR_TYPE) * (n + 1));
    memcpy(cscColIncr, cscColPtr, sizeof(MAT_PTR_TYPE) * (n + 1));

    // insert nnz to csc
    for (int row = 0; row < m; row++)
    {
        if (csrRowPtr != NULL)
        {
            for (MAT_PTR_TYPE j = csrRowPtr[row] - csrRowPtr[0]; j < csrRowPtr[row + 1] - csrRowPtr[0]; j++)
            {
                int col = csrColIdx[j] - csrRowPtr[0];

                if (cscRowIdx != NULL && cscColIncr != NULL)
                {
                    cscRowIdx[cscColIncr[col]] = row + csrRowPtr[0];
                }
                if (csrVal != NULL && cscColIncr != NULL)
                {
                    cscVal[cscColIncr[col]] = csrVal[j];
                }
                if (cscColIncr != NULL)
                {
                    cscColIncr[col]++;
                }
            }
        }
    }
    for (int i = 0; i <= n; ++i)
    {
        if (cscColPtr != NULL)
            cscColPtr[i] += csrRowPtr[0];
    }
    free(cscColIncr);
    cscColIncr = NULL;
}
double matrix_variance(double *matrix, int mA)
{
    int i, j;
    double mean = 0, variance = 0;

    // Calculate the mean of the matrix
    for (i = 0; i < mA; i++)
    {
        mean += *(matrix + i);
    }
    mean /= mA;

    // Calculate the sum of the squares of the differences between each element and the mean
    for (i = 0; i < mA; i++)
    {
        double diff = *(matrix + i) - mean;
        variance += diff * diff;
    }

    // Divide the result by the total number of elements in the matrix
    variance /= mA;

    return variance;
}
void out(char *mat_name, char *mat_group_name,
         int mA,
         int nA,
         int nnzA,
         double avg_row_A, int min_row_A, int max_row_A, double A_row_variation,
         int mB,
         int nB,
         int nnzB,
         double avg_row_B, int min_row_B, int max_row_B, double B1_row_variation,
         double avg_col_B, int min_col_B, int max_col_B, double B_row_variation, double C_row_variation,
         int nnzCub, double avg_nnzCub, int max_nnzCub, int min_nnzCub,
         double avg_time_esc,
         double min_time_esc,
         double avg_gflops_esc,
         double max_gflops_esc,
         double avg_time_hash,
         double min_time_hash,
         double avg_gflops_hash,
         double max_gflops_hash,
         double avg_time_spa,
         double min_time_spa,
         double avg_gflops_spa,
         double max_gflops_spa,
         int repeat_value)

{
    char perfname[200] = "";
    FILE *fout;
    // 记录矩阵基本信息、运算时间、性能等数据
    // 输出到文件
    strcat(perfname, "./");
    strcat(perfname, "finalresult.csv");
    fout = fopen(perfname, "a");
    fputs(mat_group_name, fout);
    fprintf(fout, ",");
    fputs(mat_name, fout);
    fprintf(fout, "%d", repeat_value);
    fprintf(fout, ",%d,%d,%d,%.4lf,%d,%d,%lf,", mA, nA, nnzA, avg_row_A, min_row_A, max_row_A, A_row_variation);
    fprintf(fout, "%d,%d,%d,%.4lf,%d,%d,%lf,", mA, nA, nnzA, avg_row_A, min_row_A, max_row_A, B1_row_variation);
    fprintf(fout, "%.4lf,%d,%d,%lf,", avg_col_B, min_col_B, max_col_B, B_row_variation);
    fprintf(fout, "%d,%.4lf,%d,%d,%lf,", nnzCub, avg_nnzCub, min_nnzCub, max_nnzCub, C_row_variation);
    fprintf(fout, "%.4f,%.4f,%.4f,%.4f,", avg_time_esc, min_time_esc, avg_gflops_esc, max_gflops_esc);
    fprintf(fout, "%.4f,%.4f,%.4f,%.4f,", avg_time_hash, min_time_hash, avg_gflops_hash, max_gflops_hash);
    fprintf(fout, "%.4f,%.4f,%.4f,%.4f", avg_time_spa, min_time_spa, avg_gflops_spa, max_gflops_spa);
    fprintf(fout, "%c", '\n');
    // 输出到屏幕
    printf("**************************************************************\n");
    printf("The Input Matrix is (%d,%d), NNZ=%d\n", mA, nA, nnzA);
    printf("%-30s\t%-15s\n", "Matrix group Name:", mat_group_name);
    printf("%-30s\t%-15s\n", "Matrix Name:", mat_name);
    printf("Test repeated: %d times\nrowA:%d\ncolA:%d\nnnzA:%d\navg_row_A:%.4lf\nmin_row_A:%d\nmax_row_A:%d\nA_row_variation:%lf\n",
           repeat_value, mA, nA, nnzA, avg_row_A, min_row_A, max_row_A, A_row_variation);

    printf("rowB:%d\ncolB:%d\nnnzB:%d\navg_row_B:%.4lf\nmin_row_B:%d\nmax_row_B:%d\nB_row_variation:%lf\navg_col_B:%.4lf\nmin_col_B:%d\nmax_col_B:%d\nB_col_variation:%lf\n",
           mA, nA, nnzA, avg_row_A, min_row_A, max_row_A, B1_row_variation, avg_col_B, min_col_B, max_col_B, B_row_variation);

    printf("nnzCub:%d\navg_nnzCub:%.4lf\nmin_nnzCub:%d\nmax_nnzCub:%d\nCub_row_variation:%lf\n", nnzCub, avg_nnzCub, min_nnzCub, max_nnzCub, C_row_variation);

    printf("avg_time_esc:%.4lf\nmin_time_esc:%.4lf\navg_gflops_esc:%.4lf\nmax_gflops_esc:%.4lf\n", avg_time_esc, min_time_esc, avg_gflops_esc, max_gflops_esc);
    printf("avg_time_hash:%.4lf\nmin_time_hash:%.4lf\navg_gflops_hash:%.4lf\nmax_gflops_hash:%.4lf\n", avg_time_hash, min_time_hash, avg_gflops_hash, max_gflops_hash);
    printf("avg_time_spa:%.4lf\nmin_time_spa:%.4lf\navg_gflops_spa:%.4lf\nmax_gflops_spa:%.4lf\n", avg_time_spa, min_time_spa, avg_gflops_spa, max_gflops_spa);
    printf("**************************************************************\n\n");
    fclose(fout);
}
#endif
