#ifndef _1_CHECK_
#define _1_CHECK_
#include "Check_d_utils.h"
#include <string.h>

int check_d_csrmultcsr(int option,
                       int nA,
                       int mA,
                       int nB,
                       int mB,
                       int *csrRowPtrA_1based,
                       int *csrColIdxA_1based,
                       double *csrValA,
                       int *csrRowPtrB_1based,
                       int *csrColIdxB_1based,
                       double *csrValB,
                       int *csrRowPtrC_1based,
                       int *csrColIdxC_1based,
                       double *csrValC,
                       int nnz)
{
    int res = 99;
    if (option == 1 || option == 0)
    {
        int request;
        int nnzC=0;
        int *rowPtrC;
        int *colIdxC;
        double *valueC;
        rowPtrC = (int *)malloc((mA + 1) * sizeof(int));
        memset(rowPtrC, 0, (mA + 1) * sizeof(int));
        request = 1;


        d_hash_csrmultcsr(&request, &mA, &nA, &nB, 
                          csrValA,csrColIdxA_1based, csrRowPtrA_1based, 
                          csrValB,csrColIdxB_1based, csrRowPtrB_1based, 
                          valueC,colIdxC, rowPtrC);
        nnzC = rowPtrC[mA] - rowPtrC[0];
        colIdxC = (int *)malloc(nnzC * sizeof(int));
        memset(colIdxC, 0, nnzC * sizeof(int));
        valueC = (double *)malloc(nnzC * sizeof(double));
        memset(valueC, 0, nnzC * sizeof(double));
        request = 2;
        d_hash_csrmultcsr(&request, &mA, &nA, &nB, 
                          csrValA,csrColIdxA_1based, csrRowPtrA_1based,
                          csrValB,csrColIdxB_1based, csrRowPtrB_1based,
                          valueC,colIdxC, rowPtrC);
        res = d_check_csr(mA, nB, 
                          nnz, csrRowPtrC_1based, csrColIdxC_1based, csrValC,
                          nnzC, rowPtrC, colIdxC, valueC);

        free(rowPtrC);
        free(colIdxC);
        free(valueC);

    }
    
    return res;
}
#endif
