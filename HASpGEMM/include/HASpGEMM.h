#ifndef _HASPGEMM_
#define _HASPGEMM_
#include "until.h"

#define setbit(x, y) x |= (1 << y)  // set the yth bit of x is 1
#define clrbit(x, y) x &= ~(1 << y) // set the yth bit of x is 0
#define getbit(x, y) ((x) >> (y)&1) // get the yth bit of x

#define Bound 256
#define binbd1 100
#define binbd2 5000

void haspgemm(const int *request,
                const int *m, const int *n, const int *k,
                double *a, int *ja, int *ia,
                double *b, int *jb, int *ib,
                double *c, int *jc, int *ic,
                const int *nzmax,int *thread_row,int core_num,int *iCub)
{
    if(*m <= Bound){ //SPA
        // step 1: compute and get the number of nnzC
        if (*request == 1)
        {
    #pragma omp parallel for
            for(int id = 0;id < core_num;id++)
            {
                for(int i = thread_row[id];i < thread_row[id+1];i++){
                    int pid = 0;
                    if (*m > 100000)
                        pid = i + 1;
                    else
                        pid = i;
                    int rowid = i;
                    int rowsize = iCub[rowid];
                    if (rowsize == 0)
                        continue;
                    int len = (*k + 31) / 32;
                    unsigned int *mask = (unsigned int *)malloc(sizeof(unsigned int) * len);
                    memset(mask, 0, sizeof(unsigned int) * len);

                    for (int offsetA = ia[rowid] - ia[0]; offsetA < ia[rowid + 1] - ia[0]; offsetA++)
                    {
                        int col = ja[offsetA] - ia[0];
                        for (int l = ib[col] - ib[0]; l < ib[col + 1] - ib[0]; l++)
                        {
                            const int key = jb[l] - ib[0];
                            setbit(mask[key / 32], key % 32);
                        }
                    }
                    int nnzr = 0;
                    for (int cid = 0; cid < *k; cid++)
                    {
                        if (getbit(mask[cid / 32], cid % 32) == 1)
                        {
                            nnzr++;
                        }
                    }
                    ic[pid] = nnzr;
                    free(mask);
                }
            }
            if (*m > 100000)
                scan_par(ic, *m + 1);
            else
                exclusive_scan(ic, *m + 1);
            for (int i = 0; i < *m + 1; i++)
            {
                ic[i] += ia[0];
            }
        }

        // step2: compute the value of C
        else
        {
            #pragma omp parallel for
            for(int id = 0;id < core_num;id++)
            {

                for(int i = thread_row[id];i < thread_row[id+1];i++){
                    int rowid = i;
                    int rowsize = iCub[rowid];
                    if (rowsize == 0)
                        continue;
                    if (rowsize == 0)
                        continue;
                    int len = (*k + 31) / 32;
                    unsigned int *mask = (unsigned int *)malloc(sizeof(unsigned int) * len);
                    memset(mask, 0, sizeof(unsigned int) * len);
                    double *d_dense_row_value = (double *)malloc((*k) * sizeof(double));
                    memset(d_dense_row_value, 0, (*k) * sizeof(double));

                    for (int offsetA = ia[rowid] - ia[0]; offsetA < ia[rowid + 1] - ia[0]; offsetA++)
                    {
                        int col = ja[offsetA] - ia[0];
                        for (int l = ib[col] - ib[0]; l < ib[col + 1] - ib[0]; l++)
                        {
                            const int key = jb[l] - ib[0];
                            setbit(mask[key / 32], key % 32);
                            d_dense_row_value[key] += b[l] * a[offsetA];
                        }
                    }

                    int nnzr = ic[rowid] - ic[0];
                    for (int cid = 0; cid < *k; cid++)
                    {
                        if (getbit(mask[cid / 32], cid % 32) == 1)
                        {
                            c[nnzr] = d_dense_row_value[cid];
                            jc[nnzr ++] = cid + ic[0];
                        }
                    }

                    free(mask);
                    free(d_dense_row_value);
                }
            }
        }
    }
    
    else{//use bin methed
        if(*request == 1){
            #pragma omp parallel for
            for(int id = 0;id < core_num;id++){
                for (int i = thread_row[id]; i < thread_row[id + 1]; i++)
                {
                    int pid = 0;
                    if (*m > 100000)
                        pid = i + 1;
                    else
                        pid = i;
                    int rowid = i;
                    int rowsize = iCub[rowid];
                    if (rowsize == 0)
                        continue;
                    
                    if(rowsize <= binbd1){
                        //esc
                        int *jCub = (int *)malloc(rowsize * sizeof(int));
                        memset(jCub, 0, rowsize * sizeof(int));

                        int incr = 0;
                        for (int l = ia[rowid] - ia[0]; l < ia[rowid + 1] - ia[0]; l++)
                        {
                            int rowB = ja[l] - ia[0];
                            for (int k = ib[rowB] - ib[0]; k < ib[rowB + 1] - ib[0]; k++)
                            {
                                jCub[incr] = jb[k] - ib[0];
                                incr++;
                            }
                        }
                        quick_sort_key1(&jCub[0], rowsize);

                        int nnzr = rowsize > 0 ? 1 : 0;
                        for (int idx = 1; idx < rowsize; idx++)
                        {
                            nnzr = jCub[idx] == jCub[idx - 1] ? nnzr : nnzr + 1;
                        }
                        ic[pid] = nnzr;
                        free(jCub);
                    }
                    else if(rowsize <= binbd2){
                        //hash
                        int hashsize_full_reg = rowsize / 0.75;
                        int *tmpHashtable = (int *)malloc(hashsize_full_reg * sizeof(int));
                        memset(tmpHashtable, -1, sizeof(int) * hashsize_full_reg);
                        for (int blkj = ia[rowid] - ia[0]; blkj < ia[rowid + 1] - ia[0]; blkj++)
                        {
                            int col = ja[blkj] - ia[0];
                            for (int l = ib[col] - ib[0]; l < ib[col + 1] - ib[0]; l++)
                            {
                                const int key = jb[l] - ib[0];
                                int hashadr = (key) % hashsize_full_reg;
                                while (1)
                                {
                                    const int keyexist = tmpHashtable[hashadr];
                                    if (keyexist == key)
                                    {
                                        break;
                                    }
                                    else if (keyexist == -1)
                                    {
                                        tmpHashtable[hashadr] = key;
                                        ic[pid]++;
                                        break;
                                    }
                                    else
                                    {
                                        hashadr = (hashadr + 1) % hashsize_full_reg;
                                    }
                                }
                            }
                        }
                        free(tmpHashtable);
                    }
                    else{
                        //spa
                        int len = (*k + 31) / 32;
                        unsigned int *mask = (unsigned int *)malloc(sizeof(unsigned int) * len);
                        memset(mask, 0, sizeof(unsigned int) * len);
                        for (int offsetA = ia[rowid] - ia[0]; offsetA < ia[rowid + 1] - ia[0]; offsetA++)
                        {
                            int col = ja[offsetA] - ia[0];
                            for (int l = ib[col] - ib[0]; l < ib[col + 1] - ib[0]; l++)
                            {
                                const int key = jb[l] - ib[0];
                                setbit(mask[key / 32], key % 32);
                            }
                        }
                        int nnzr = 0;
                        for (int cid = 0; cid < *k; cid++)
                        {
                            if (getbit(mask[cid / 32], cid % 32) == 1)
                            {
                                nnzr++;
                            }
                        }
                        ic[pid] = nnzr;
                        free(mask);
                    }
                }
            }
            if (*m > 100000)
                scan_par(ic, *m +1);
            else
                exclusive_scan(ic, *m +1);
            for (int i = 0; i < *m+1; i++)
            {
                ic[i] += ia[0];
            }
        }
        else{
            #pragma omp parallel for
            for(int id = 0;id < core_num;id++){
                for (int i = thread_row[id]; i < thread_row[id + 1]; i++)
                {
                    int rowid = i;
                    int rowsize = iCub[rowid];
                    if (rowsize == 0)
                        continue;

                    if(rowsize <= binbd1){
                        //esc
                        int *jCub = (int *)malloc(rowsize * sizeof(int));
                        double *valCub = (double *)malloc(rowsize * sizeof(double));
                        int *d_flagCub = (int *)malloc(rowsize * sizeof(int));
                        memset(jCub, 0, rowsize * sizeof(int));
                        memset(valCub, 0, rowsize * sizeof(double));
                        memset(d_flagCub, 0, rowsize * sizeof(int));

                        int incr = 0;
                        for (int l = ia[rowid] - ia[0]; l < ia[rowid + 1] - ia[0]; l++)
                        {
                            int rowB = ja[l] - ia[0];
                            double val = a[l];
                            for (int k = ib[rowB] - ib[0]; k < ib[rowB + 1] - ib[0]; k++)
                            {
                                jCub[incr] = jb[k];
                                valCub[incr] = val * b[k];
                                incr++;
                            }
                        }

                        // sort
                        quick_sort_key_val_pair1(jCub, valCub, rowsize);

                        // compress
                        d_flagCub[0] = 1;
                        for (int idx = 0; idx < rowsize - 1; idx++)
                        {
                            d_flagCub[idx + 1] = jCub[idx + 1] == jCub[idx] ? 0 : 1;
                        }
                        segmented_sum(valCub, d_flagCub, rowsize);

                        int incrn = 0;
                        for (int idx = 0; idx < rowsize; idx++)
                        {
                            if (d_flagCub[idx] == 1)
                            {
                                jc[ic[rowid] - ic[0] + incrn] = jCub[idx];
                                c[ic[rowid] - ic[0] + incrn] = valCub[idx];
                                incrn++;
                            }
                        }

                        free(jCub);
                        free(valCub);
                        free(d_flagCub);
                    }
                    else if(rowsize <= binbd2){
                        //hash
                        int hashsize_full_reg = (ic[rowid + 1] - ic[rowid]) / 0.75;
                        int *tmpHashtable = (int *)malloc(hashsize_full_reg * sizeof(int));
                        memset(tmpHashtable, -1, sizeof(int) * hashsize_full_reg);
                        double *tmpValue = (double *)malloc(hashsize_full_reg * sizeof(double));
                        memset(tmpValue, 0, sizeof(double) * hashsize_full_reg);
                        for (int blkj = ia[rowid] - ia[0]; blkj < ia[rowid + 1] - ia[0]; blkj++)
                        {
                            int col = ja[blkj] - ia[0];
                            for (int l = ib[col] - ib[0]; l < ib[col + 1] - ib[0]; l++)
                            {
                                const int key = jb[l] - ib[0];
                                int hashadr = (key * 107) % hashsize_full_reg;
                                while (1)
                                {
                                    const int keyexist = tmpHashtable[hashadr];
                                    if (keyexist == key)
                                    {
                                        tmpValue[hashadr] += b[l] * a[blkj];
                                        break;
                                    }
                                    else if (keyexist == -1)
                                    {
                                        tmpHashtable[hashadr] = key;
                                        tmpValue[hashadr] = b[l] * a[blkj];
                                        break;
                                    }
                                    else
                                    {
                                        hashadr = (hashadr + 1) % hashsize_full_reg;
                                    }
                                }
                            }
                        }
                        int cptr = ic[rowid] - ic[0];
                        for (int k = 0; k < hashsize_full_reg; k++)
                        {
                            if (tmpHashtable[k] != -1)
                            {
                                jc[cptr] = tmpHashtable[k] + ib[0];
                                c[cptr] = tmpValue[k];
                                cptr++;
                            }
                        }
                        free(tmpHashtable);
                        free(tmpValue);
                        int nnzcnt = ic[rowid + 1] - ic[rowid];
                        quick_sort_key_val_pair1(jc + ic[rowid] - ic[0], c + ic[rowid] - ic[0], nnzcnt);
                    }
                    else{
                        //spa
                        int len = (*k + 31) / 32;
                        unsigned int *mask = (unsigned int *)malloc(sizeof(unsigned int) * len);
                        memset(mask, 0, sizeof(unsigned int) * len);
                        double *d_dense_row_value = (double *)malloc((*k) * sizeof(double));
                        memset(d_dense_row_value, 0, (*k) * sizeof(double));

                        for (int offsetA = ia[rowid] - ia[0]; offsetA < ia[rowid + 1] - ia[0]; offsetA++)
                        {
                            int col = ja[offsetA] - ia[0];
                            for (int l = ib[col] - ib[0]; l < ib[col + 1] - ib[0]; l++)
                            {
                                const int key = jb[l] - ib[0];
                                setbit(mask[key / 32], key % 32);
                                d_dense_row_value[key] += b[l] * a[offsetA];
                            }
                        }

                        int nnzr = ic[rowid] - ic[0];
                        for (int cid = 0; cid < *k; cid++)
                        {
                            if (getbit(mask[cid / 32], cid % 32) == 1)
                            {
                                c[nnzr] = d_dense_row_value[cid];
                                jc[nnzr] = cid + ic[0];
                                nnzr++;
                            }
                        }

                        free(mask);
                        free(d_dense_row_value);
                    }
                }
            }
            
        }
    }
}
#endif
