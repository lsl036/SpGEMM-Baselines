#include <omp.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#ifndef MAT_PTR_TYPE
#define MAT_PTR_TYPE int
#endif
#ifndef MAT_VAL_TYPE
#define MAT_VAL_TYPE double
#endif

#ifndef P_NUM
#define P_NUM 16
#endif
#ifndef E_NUM
#define E_NUM 8
#endif