#include <stdio.h>
#include <stdlib.h>
void d_swap_key(int *a, int *b)
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

void d_swap_val(double *a, double *b)
{
    double tmp = *a;
    if (a != NULL && b != NULL)
    {
        *a = *b;
    }
    if (b != NULL)
    {
        *b = tmp;
    }
}

int d_partition_key_val_pair(int *key, double *val, int length, int pivot_index)
{
    int i = 0;
    int small_length = pivot_index;

    int pivot = key[pivot_index];
    if (val != NULL && key != NULL)
    {
        d_swap_key(&key[pivot_index], &key[pivot_index + (length - 1)]);
        d_swap_val(&val[pivot_index], &val[pivot_index + (length - 1)]);
    }

    for (; i < length; i++)
    {
        if (key != NULL)
        {
            if (key[pivot_index + i] < pivot)
            {
                d_swap_key(&key[pivot_index + i], &key[small_length]);
                if (val != NULL)
                {
                    d_swap_val(&val[pivot_index + i], &val[small_length]);
                }
                small_length++;
            }
        }
    }

    if (key != NULL)
    {
        d_swap_key(&key[pivot_index + length - 1], &key[small_length]);
    }
    if (val != NULL)
    {
        d_swap_val(&val[pivot_index + length - 1], &val[small_length]);
    }

    return small_length;
}

void d_quick_sort_key_val_pair(int *key, double *val, int length)
{
    if (length == 0 || length == 1)
    {
        return;
    }

    int small_length = d_partition_key_val_pair(key, val, length, 0);
    d_quick_sort_key_val_pair(key, val, small_length);
    if (val != NULL && key != NULL)
    {
        d_quick_sort_key_val_pair(&key[small_length + 1], &val[small_length + 1], length - small_length - 1);
    }
}

void d_exclusive_scan(int *input, int length)
{
    if (length == 0 || length == 1)
        return;

    int old_val, new_val;

    old_val = input[0];
    input[0] = 0;
    for (int i = 1; i < length; i++)
    {
        new_val = input[i];
        input[i] = old_val + input[i - 1];
        old_val = new_val;
    }
}

void d_hash_csrmultcsr(const int *request,
                       const int *m, const int *n, const int *k,
                       double *a, int *ja, int *ia,
                       double *b, int *jb, int *ib,
                       double *c, int *jc, int *ic)
{
    int nthreads = omp_get_max_threads();
    int hashsize_full_reg = 0;
    for (int blki = 0; blki < *m; blki++)
    {
        int max = 0;
        for (int l = ia[blki] - ia[0]; l < ia[blki + 1] - ia[0]; l++)
        {
            int cola = ja[l] - ia[0];
            max += ib[cola + 1] - ib[cola];
        }
        if (max > hashsize_full_reg)
            hashsize_full_reg = max;
    }
    int *tmpIdx2D0_g = (int *)malloc(nthreads * hashsize_full_reg * sizeof(int));
    if (*request == 1)
    {
#pragma omp parallel for
        for (int blki = 0; blki < *m; blki++)
        {
            int thread_id = omp_get_thread_num();

            int *tmpIdx2D0 = tmpIdx2D0_g + thread_id * hashsize_full_reg;
            for (int l = 0; l < hashsize_full_reg; l++) // hashsize_full_reg is the hash length assigned to the thread
            {
                tmpIdx2D0[l] = -1;
            } // 0x80000000=-1
            if (tmpIdx2D0 != NULL)
            {
                for (int blkj = ia[blki] - ia[0]; blkj < ia[blki + 1] - ia[0]; blkj++)
                {
                    int col = ja[blkj] - ia[0];
                    for (int l = ib[col] - ib[0]; l < ib[col + 1] - ib[0]; l++)
                    {
                        const int key = jb[l];
                        int hashadr = (key) % hashsize_full_reg;
                        while (1)
                        {
                            const int keyexist = tmpIdx2D0[hashadr]; // tmpIdx2Dthread[hashadr];
                            if (keyexist == key)
                            {
                                break;
                            }
                            else if (keyexist == -1)
                            {
                                tmpIdx2D0[hashadr] = key;
                                ic[blki]++;
                                break;
                            }
                            else
                            {
                                hashadr = (hashadr + 1) % hashsize_full_reg;
                                // in step 1, it is not possible to overflow, since the assigned space is upper bound
                            }
                        }
                    }
                }
            }
            for (int l = 0; l < hashsize_full_reg; l++) // hashsize_full_reg is the hash length assigned to the thread
                tmpIdx2D0[l] = -1;                      // 0x80000000=-1 ,pos is the starting position of the thread
        }
        d_exclusive_scan(ic, *m + 1);
        for (int i = 0; i < *m + 1; i++)
        {
            ic[i] += ia[0];
        }
    }
    else if (*request == 2)
    {
        MAT_VAL_TYPE *tmpVal2D0_g = (MAT_VAL_TYPE *)malloc(hashsize_full_reg * nthreads * sizeof(MAT_VAL_TYPE));
        memset(tmpVal2D0_g, 0, hashsize_full_reg * nthreads * sizeof(MAT_VAL_TYPE));
#pragma omp parallel for
        for (int blki = 0; blki < *m; blki++)
        {
            int thread_id = omp_get_thread_num();

            int *tmpIdx2D0 = tmpIdx2D0_g + thread_id * hashsize_full_reg;
            MAT_VAL_TYPE *tmpVal2D0 = tmpVal2D0_g + thread_id * hashsize_full_reg;
            for (int l = 0; l < hashsize_full_reg; l++) // hashsize_full_reg is the hash length assigned to the thread
                tmpIdx2D0[l] = -1;                      // 0x80000000=-1

            for (int blkj = ia[blki] - ia[0]; blkj < ia[blki + 1] - ia[0]; blkj++)
            {
                int col = ja[blkj] - ia[0];
                for (int l = ib[col] - ib[0]; l < ib[col + 1] - ib[0]; l++)
                {
                    const int key = jb[l];
                    int hashadr = (key * 107) % hashsize_full_reg;
                    while (1)
                    {
                        const int keyexist = tmpIdx2D0[hashadr]; // tmpIdx2Dthread[hashadr];
                        if (keyexist == key)
                        {
                            tmpVal2D0[hashadr] += b[l] * a[blkj];
                            break;
                        }
                        else if (keyexist == -1)
                        {
                            tmpIdx2D0[hashadr] = key;
                            tmpVal2D0[hashadr] = b[l] * a[blkj];
                            break;
                        }
                        else
                        {
                            hashadr = (hashadr + 1) % hashsize_full_reg;
                            // in step 1, it is not possible to overflow, since the assigned space is upper bound
                        }
                    }
                }
            }
            int cptr = ic[blki] - ic[0];
            for (int k = 0; k < hashsize_full_reg; k++)
            {
                if (tmpIdx2D0[k] != -1)
                {
                    jc[cptr] = tmpIdx2D0[k];
                    c[cptr] = tmpVal2D0[k];
                    cptr++;
                }
            }

            for (int l = 0; l < hashsize_full_reg; l++) // hashsize_full_reg is the hash length assigned to the thread
                tmpIdx2D0[l] = -1;                      // 0x80000000=-1 ,pos is the starting position of the thread
            memset(tmpVal2D0, 0, hashsize_full_reg * sizeof(MAT_VAL_TYPE));
        }

        free(tmpIdx2D0_g);

        free(tmpVal2D0_g);

#pragma omp parallel for
        for (int i = 0; i < *m; i++)
        {
            int nnzcnt = ic[i + 1] - ic[i];
            d_quick_sort_key_val_pair(jc + ic[i] - ic[0], c + ic[i] - ic[0], nnzcnt);
        }
    }
}

void d_csrmultd(const int *m, const int *n, const int *k,
                double *a, int *ja, int *ia,
                double *b, int *jb, int *ib, double *c)
{
    int nthreads = omp_get_max_threads();
    nthreads = nthreads > *m ? *m : nthreads;
    int threadrow = *m / nthreads;

#pragma omp parallel for
    for (int tid = 0; tid < nthreads; tid++)
    {
        int hashsize_full_reg = 0;
        int rowstart = tid * threadrow;
        int rowend;
        if (tid < nthreads - 1)
        {
            rowend = (tid + 1) * threadrow;
        }
        else if (tid == nthreads - 1)
        {
            rowend = *m;
        }
        else
        {
        }
        for (int j = rowstart; j < rowend; j++)
        {
            for (int i = ia[j] - ia[0]; i < ia[j + 1] - ia[0]; i++)
            {
                int col = ja[i] - ia[0];
                for (int l1 = ib[col] - ib[0]; l1 < ib[col + 1] - ib[0]; l1++)
                {
                    int key = jb[l1] - ib[0];
                    c[j * *k + key] += a[i] * b[l1];
                }
            }
        }
    }
}

int d_check_dense(int m, int n, double *valA, double *valB)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (fabs(valA[i * n + j] - valB[i * n + j]) > 1e-6)
            {
                printf("i=%d,j=%d,A=%lf,B=%lf,the value is dfferent, may some error!\n", i, j, valA[i * n + j], valB[i * n + j]);
                return 1;
            }
        }
    }

    return 0;
}

int d_check_csr(int m, int n,
                int nnzA, int *rptA, int *cidA, double *valA,
                int nnzB, int *rptB, int *cidB, double *valB)
{
    if (nnzA != nnzB)
        printf("The nnzc is different! may have some error!\n");

    for (int i = 0; i < m + 1; i++)
    {
        if (rptA[i] != rptB[i])
        {
            printf("The rowptr is different! may some error in line %d!\n", i);
            return 1;
        }
    }

    for (int i = 0; i < nnzA; i++)
    {
        // printf("%d %d\n",cidA[i],cidB[i]);
        if (cidA[i] != cidB[i])
        {
            printf("The colid is different! may some error ! i=%d, cidA[i] = %d, cidB[i] = %d\n", i, cidA[i], cidB[i]);
            return 2;
        }
        if (fabs(valA[i] - valB[i]) > 1e-6)
        {
            printf("The value is different! may some error !i=%d, valA[i] = %lf, valB[i] = %lf\n", i, valA[i], valB[i]);
            return 3;
        }
    }

    return 0;
}

void d_clear_upper(double *B, int m)
{
    for (int i = 1; i < m; ++i)
    {
        memset(B + i * m, 0, (i) * sizeof(double));
    }
}

void d_csrCpToUpTri(int m, int n,
                    const int *csrRowPtr, const int *csrColIdx, const double *csrVal,
                    int **uptRowPtr1, int **uptColIdx1, double **uptVal1)
{
    int *uptRowPtr = (int*)malloc((m + 1) * sizeof(int));
    memset(uptRowPtr, 0, (m + 1) * sizeof(int));
    int nnz = csrRowPtr[m] - csrRowPtr[0];
    int *uptColIdx = (int*)malloc((nnz + m) * sizeof(int));
    double *uptVal = (double*)malloc((nnz + m) * sizeof(double));

    int nnz_upt = 0;
    uptRowPtr[0] = csrRowPtr[0];
    for (int i = 0; i < m; ++i)
    {
        if (csrRowPtr != NULL)
            for (int j = csrRowPtr[i] - csrRowPtr[0]; j < csrRowPtr[i + 1] - csrRowPtr[0]; ++j)
            {
                if (csrColIdx != NULL)
                {
                    if (csrColIdx[j] - csrRowPtr[0] < i)
                        continue;
                    else
                    {
                        uptVal[nnz_upt] = csrVal[j];
                        uptColIdx[nnz_upt] = csrColIdx[j];
                        ++nnz_upt;
                    }
                }
            }
        uptRowPtr[i + 1] = nnz_upt + csrRowPtr[0];
    }
    *uptRowPtr1 = uptRowPtr;
    *uptColIdx1 = uptColIdx;
    *uptVal1 = uptVal;
}

void d_matrix_transposition(const int m,
                            const int n,
                            const int nnz,
                            const int *csrRowPtr,
                            const int *csrColIdx,
                            const double *csrVal,
                            int *cscRowIdx,
                            int *cscColPtr,
                            double *cscVal)
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
    d_exclusive_scan(cscColPtr, n + 1);

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

void d_csrLowToFull(int m, int n,
                    const int *lowRowPtr, const int *lowColIdx, const double *lowVal,
                    int **fullRowPtr1, int **fullColIdx1, double **fullVal1)
{
    int low_nnz = lowRowPtr[m] - lowRowPtr[0];
    int *upRowPtr = (int *)malloc((m + 1) * sizeof(int));
    int *upColIdx = (int *)malloc(low_nnz * sizeof(int));
    double *upVal = (double *)malloc(low_nnz * sizeof(double));
    memset(upRowPtr, 0, (m + 1) * sizeof(int));
    memset(upColIdx, 0, low_nnz * sizeof(int));
    memset(upVal, 0, low_nnz * sizeof(double));

    d_matrix_transposition(m, n, low_nnz,
                           lowRowPtr, lowColIdx, lowVal,
                           upColIdx, upRowPtr, upVal);

    int *fullRowPtr = (int *)malloc((m + 1) * sizeof(int));
    int *fullColIdx = (int *)malloc(low_nnz * 2 * sizeof(int));
    double *fullVal = (double *)malloc(low_nnz * 2 * sizeof(double));
    memset(fullRowPtr, 0, (m + 1) * sizeof(int));
    memset(fullColIdx, 0, low_nnz * 2 * sizeof(int));
    memset(fullVal, 0, low_nnz * 2 * sizeof(double));

    int nnz_full = 0;
    fullRowPtr[0] = lowRowPtr[0];
    for (int i = 0; i < m; i++)
    {
        int cur_nnz = lowRowPtr[i + 1] - lowRowPtr[i];
        int begJ = lowRowPtr[i] - lowRowPtr[0];
        fullRowPtr[i + 1] += cur_nnz;
        memcpy(fullColIdx + nnz_full, lowColIdx + begJ, sizeof(int) * (cur_nnz));
        memcpy(fullVal + nnz_full, lowVal + begJ, sizeof(double) * (cur_nnz));
        nnz_full += cur_nnz;

        if (i + upRowPtr[0] == upColIdx[upRowPtr[i] - upRowPtr[0]])
        {
            cur_nnz = upRowPtr[i + 1] - upRowPtr[i] - 1;
            begJ = upRowPtr[i] - upRowPtr[0] + 1;
        }
        else
        {
            cur_nnz = upRowPtr[i + 1] - upRowPtr[i];
            begJ = upRowPtr[i] - upRowPtr[0];
        }

        fullRowPtr[i + 1] += cur_nnz;
        memcpy(fullColIdx + nnz_full, upColIdx + begJ, sizeof(int) * (cur_nnz));
        memcpy(fullVal + nnz_full, upVal + begJ, sizeof(double) * (cur_nnz));
        nnz_full += cur_nnz;
        fullRowPtr[i + 1] = nnz_full + fullRowPtr[0];
    }
    *fullRowPtr1 = fullRowPtr;
    *fullColIdx1 = fullColIdx;
    *fullVal1 = fullVal;

    free(upRowPtr);
    free(upColIdx);
    free(upVal);
}

void d_csr_mult_d(int m, int n, int k,
                  const int *ia, const int *ja, const double *val,
                  const double *B, double *res, double alpha)
{
#pragma omp parallel for
    for (int i = 0; i < m; ++i)
    {
        for (int kk = 0; kk < k; ++kk)
        {
            double sum = 0;
            for (int j = ia[i] - ia[0]; j < ia[i + 1] - ia[0]; ++j)
            {
                sum += val[j] * B[(ja[j] - ia[0]) * k + kk];
            }
            res[i * k + kk] += alpha * sum;
        }
    }
}

void d_d_mult_csrT(int m, int n, int k,
                   const double *B,
                   const int *tia, const int *tja, const double *tval,
                   double *res)
{
#pragma omp parallel for
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < k; j++)
        {
            double val = B[i * k + j];
            int colidx = j;
            for (int kk = tia[colidx] - tia[0]; kk < tia[colidx + 1] - tia[0]; kk++)
            {
                int colb = tja[kk];
                double valc = tval[kk];
                res[i * m + colb] += val * valc;
            }
        }
    }
}
