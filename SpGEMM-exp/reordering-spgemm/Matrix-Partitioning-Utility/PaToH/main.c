#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <stdbool.h>
#include <sys/stat.h>

#include "patoh.h"
#include "smsh.c"
#include "mmio.c"

void PrintHelp(const char *program_name)
{
    printf("PaToH Matrix Partitioning Tool\n");
    printf("================================\n\n");
    printf("Usage: %s <matrix-file> <partition-type> <num-partitions> <quality> <seed> <save-to-file>\n\n", program_name);
    printf("Arguments:\n");
    printf("  <matrix-file>          Path to the Matrix Market format file (.mtx)\n");
    printf("  <partition-type>       Partitioning objective:\n");
    printf("                           - cutpart  : Minimize cutsize (edge-cut minimization)\n");
    printf("                           - conpart  : Minimize connectivity (connectivity minimization)\n");
    printf("  <num-partitions>       Number of partitions (k-way partitioning, e.g., 64)\n");
    printf("  <quality>              Partitioning quality/speed trade-off:\n");
    printf("                           - quality  : High quality (slower, better results)\n");
    printf("                           - default  : Balanced quality and speed\n");
    printf("                           - speed    : Fast partitioning (faster, lower quality)\n");
    printf("  <seed>                 Random seed for reproducibility (integer, e.g., 1)\n");
    printf("  <save-to-file>         Save partition vector to file:\n");
    printf("                           - 0        : Do not save partvec.txt\n");
    printf("                           - 1        : Save partvec.txt\n\n");
    printf("Output Files:\n");
    printf("  The program always generates the following files in the 'output/' directory:\n");
    printf("    - PaToH_<matrix>_<type>_k<num>_<quality>_s<seed>_partinfo.txt  : Partition information\n");
    printf("    - PaToH_<matrix>_<type>_k<num>_<quality>_s<seed>_timeinfo.txt   : Timing information\n");
    printf("  If save-to-file=1, also generates:\n");
    printf("    - PaToH_<matrix>_<type>_k<num>_<quality>_s<seed>_partvec.txt    : Partition vector\n\n");
    printf("Examples:\n");
    printf("  %s matrix.mtx cutpart 64 quality 1 1\n", program_name);
    printf("  %s matrix.mtx conpart 32 default 0 0\n\n", program_name);
    printf("Options:\n");
    printf("  -h, --help             Show this help message\n\n");
}

void PrintInfo(int _k, int *partweights, int cut, int _nconst, char *base_file_name, char *output_file_name, int non_zero_count, int row_count, int col_count)
{
    double *avg, *maxi, maxall = -1.0;
    int i, j;

    FILE *info_fp;

    info_fp = fopen(output_file_name, "w");
    fprintf(info_fp, "File Name: %s\n", base_file_name);
    fprintf(info_fp, "Non-Zero Count: %d\n", non_zero_count);
    fprintf(info_fp, "Row and Column Count (M x N): %d and %d\n", row_count, col_count);
    fprintf(info_fp, "-------------------------------------------------------------------");
    fprintf(info_fp, "\n Partitioner: %s", (_nconst > 1) ? "Multi-Constraint" : "Single-Constraint");
    fprintf(info_fp, "\n %d-way cutsize = %d \n", _k, cut);

    avg = (double *)malloc(sizeof(double) * _nconst);
    maxi = (double *)malloc(sizeof(double) * _nconst);
    for (i = 0; i < _nconst; ++i)
        maxi[i] = avg[i] = 0.0;
    for (i = 0; i < _k; ++i)
        for (j = 0; j < _nconst; ++j)
            avg[j] += partweights[i * _nconst + j];
    for (i = 0; i < _nconst; ++i)
    {
        maxi[i] = 0.0;
        avg[i] /= (double)_k;
    }

    for (i = 0; i < _k; ++i)
    {
        // fprintf(info_fp, "\n %3d :", i);
        for (j = 0; j < _nconst; ++j)
        {
            double im = (double)((double)partweights[i * _nconst + j] - avg[j]) / avg[j];

            maxi[j] = (maxi[j] > im) ? maxi[j] : im;
            //     fprintf(info_fp, "%10d ", partweights[i * _nconst + j]);
        }
    }

    for (j = 0; j < _nconst; ++j)
        maxall = (maxi[j] > maxall) ? maxi[j] : maxall;
    fprintf(info_fp, "\n MaxImbals are (as %%): %.3lf", 100.0 * maxall);
    fprintf(info_fp, "\n      ");
    for (i = 0; i < _nconst; ++i)
        fprintf(info_fp, "%10.1lf \n", 100.0 * maxi[i]);

    fprintf(info_fp, "\n PartWeights are:\n");

    for (i = 0; i < _k; ++i)
    {
        fprintf(info_fp, "\n %3d :", i);
        for (j = 0; j < _nconst; ++j)
        {
            fprintf(info_fp, "%10d ", partweights[i * _nconst + j]);
        }
    }

    fprintf(info_fp, "\n");
    fclose(info_fp);
    free(maxi);
    free(avg);
}

int main(int argc, char *argv[])
{
    clock_t program_begins = clock();
    /* Declarations for reading the matrix */

    int ret_code;
    MM_typecode matcode;
    FILE *f;
    bool save_to_file = false;
    int M, N, nz;
    int i, *I_complete, *J_complete;

    /* Declarations for measuring time */
    double elapsedTime = 0.0;
    double elapsed_time_patoh_part_only = 0.0;

    /* Processing the input */

    // Check for help option
    if (argc == 2 && (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))
    {
        PrintHelp(argv[0]);
        exit(0);
    }

    if (argc != 7)
    {
        printf("Error: Incorrect number of arguments.\n\n");
        PrintHelp(argv[0]);
        exit(1);
    }

    // Validate matrix file
        if ((f = fopen(argv[1], "r")) == NULL)
        {
        printf("Error: Could not open matrix file: %s\n", argv[1]);
        printf("Please make sure the file path is correct and the file exists.\n\n");
        PrintHelp(argv[0]);
            exit(1);
        }

    // Validate partition type
        if ((strcmp(argv[2], "conpart") != 0) && (strcmp(argv[2], "cutpart") != 0))
        {
        printf("Error: Invalid partition type: %s\n", argv[2]);
        printf("Valid options are: 'conpart' or 'cutpart'\n\n");
        PrintHelp(argv[0]);
            exit(1);
        }

    // Validate partition quality
        if ((strcmp(argv[4], "speed") != 0) && (strcmp(argv[4], "default") != 0) && (strcmp(argv[4], "quality") != 0))
        {
        printf("Error: Invalid partition quality: %s\n", argv[4]);
        printf("Valid options are: 'speed', 'default', or 'quality'\n\n");
        PrintHelp(argv[0]);
            exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Could not process Matrix Market banner.\n");
        exit(1);
    }

    if ((ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz)) != 0)
    {
        printf("Could not read matrix dimensions.\n");
        exit(1);
    }

    if ((strcmp(matcode, "MCRG") == 0) || (strcmp(matcode, "MCIG") == 0) || (strcmp(matcode, "MCPG") == 0) || (strcmp(matcode, "MCCG") == 0))
    {

        I_complete = (int *)calloc(nz, sizeof(int));
        J_complete = (int *)calloc(nz, sizeof(int));

        for (i = 0; i < nz; i++)
        {
            fscanf(f, "%d %d", &I_complete[i], &J_complete[i]);
            fscanf(f, "%*[^\n]\n");
            /* adjust from 1-based to 0-based */
            I_complete[i]--;
            J_complete[i]--;
        }
    }

    /* If the matrix is symmetric, we need to construct the other half */

    else if ((strcmp(matcode, "MCRS") == 0) || (strcmp(matcode, "MCIS") == 0) || (strcmp(matcode, "MCPS") == 0) || (strcmp(matcode, "MCCS") == 0) || (strcmp(matcode, "MCCH") == 0) || (strcmp(matcode, "MCRK") == 0))
    {
        // Use long long to avoid integer overflow when nz is large
        // In worst case, we need nz (original) + nz (symmetric) = 2 * nz
        long long alloc_size = (long long)2 * (long long)nz;
        I_complete = (int *)calloc(alloc_size, sizeof(int));
        J_complete = (int *)calloc(alloc_size, sizeof(int));
        
        if (I_complete == NULL || J_complete == NULL)
        {
            printf("Error: Memory allocation failed. Matrix too large (nz = %d, trying to allocate %lld elements).\n", nz, alloc_size);
            printf("Required memory: %.2f GB\n", (double)alloc_size * sizeof(int) / (1024.0 * 1024.0 * 1024.0));
            exit(1);
        }

        int i_index = 0;

        for (i = 0; i < nz; i++)
        {
            fscanf(f, "%d %d", &I_complete[i], &J_complete[i]);
            fscanf(f, "%*[^\n]\n");

            if (I_complete[i] == J_complete[i])
            {
                /* adjust from 1-based to 0-based */
                I_complete[i]--;
                J_complete[i]--;
            }
            else
            {
                /* adjust from 1-based to 0-based */
                I_complete[i]--;
                J_complete[i]--;
                // Use long long to avoid integer overflow in array indexing
                long long idx = (long long)nz + (long long)i_index;
                if (idx >= alloc_size)
                {
                    printf("Error: Array bounds exceeded. i_index = %d, nz = %d, idx = %lld, max allowed = %lld\n", 
                           i_index, nz, idx, alloc_size - 1);
                    exit(1);
                }
                J_complete[idx] = I_complete[i];
                I_complete[idx] = J_complete[i];
                i_index++;
            }
        }
        nz += i_index;
    }
    else
    {
        printf("This matrix type is not supported: %s \n", matcode);
        exit(1);
    }

    /* We need the values to be sorted with respect to column index to convert COO to CSC */

    if (!isSorted(I_complete, J_complete, nz))
    {
        quicksort(I_complete, J_complete, nz);
    }

    /**** Create output file names ****/

    char partvecFileName[256];
    memset(partvecFileName, 0, 256 * sizeof(char));
    strcat(strcat(strcat(strcat(strcat(strcat(strcat(strcat(strcat(strcat(strcat(partvecFileName, "output/PaToH_"), basename(argv[1])), "_"), argv[2]), "_k"), argv[3]), "_"), argv[4]), "_s"), argv[5]), "_partvec.txt");

    char partinfoFileName[256];
    memset(partinfoFileName, 0, 256 * sizeof(char));
    strcat(strcat(strcat(strcat(strcat(strcat(strcat(strcat(strcat(strcat(strcat(partinfoFileName, "output/PaToH_"), basename(argv[1])), "_"), argv[2]), "_k"), argv[3]), "_"), argv[4]), "_s"), argv[5]), "_partinfo.txt");

    char timeinfoFileName[256];
    memset(timeinfoFileName, 0, 256 * sizeof(char));
    strcat(strcat(strcat(strcat(strcat(strcat(strcat(strcat(strcat(strcat(strcat(timeinfoFileName, "output/PaToH_"), basename(argv[1])), "_"), argv[2]), "_k"), argv[3]), "_"), argv[4]), "_s"), argv[5]), "_timeinfo.txt");

    /* Part that uses PaToH */

    PaToH_Parameters args;
    int _c, _n, _nconst, *cwghts, *nwghts, *xpins, *pins, *partvec, cut, *partweights;

    _c = M;
    _n = N;
    _nconst = 1;

    cwghts = (int *)malloc(_c * sizeof(int));
    for (int i = 0; i < _c; i++)
    {
        cwghts[i] = 1;
    }

    nwghts = (int *)malloc(_n * sizeof(int));
    for (int i = 0; i < _n; i++)
    {
        nwghts[i] = 1;
    }

    // nwghts = NULL;

    /* Convert COO to CSC, and write the column pointer and row indices to xpins and pins respectivelly. */

    xpins = (int *)calloc(_n + 1, sizeof(int));
    pins = (int *)calloc(nz, sizeof(int));

    for (i = 0; i < nz; i++)
    {
        pins[i] = I_complete[i];
        xpins[J_complete[i] + 1]++;
    }
    /* Prefix sum: xpins has _n+1 elements (indices 0.._n); only iterate i=0.._n-1 to avoid writing xpins[_n+1] */
    for (i = 0; i < _n; i++)
    {
        xpins[i + 1] += xpins[i];
    }

    free(I_complete);
    free(J_complete);

    /*
        FILE *xpins_fp;

        xpins_fp = fopen("xpins_file.txt", "w");

        for (i = 0; i < N + 1; i++)
        {
            fprintf(xpins_fp, "%d\n", xpins[i]);
        }
        fclose(xpins_fp);

        FILE *pins_fp;

        pins_fp = fopen("pins_file.txt", "w");

        for (i = 0; i < nz; i++)
        {
            fprintf(pins_fp, "%d\n", pins[i]);
        }
        fclose(pins_fp);
     */

    partvec = (int *)malloc(_c * sizeof(int));
    partweights = (int *)malloc(atoi(argv[3]) * _nconst * sizeof(int));

    if (!strcmp(argv[2], "conpart"))
    {

        if (!strcmp(argv[4], "quality"))
        {
            PaToH_Initialize_Parameters(&args, PATOH_CONPART, PATOH_SUGPARAM_QUALITY);
        }
        else if (!strcmp(argv[4], "speed"))
        {
            PaToH_Initialize_Parameters(&args, PATOH_CONPART, PATOH_SUGPARAM_SPEED);
        }
        else if (!strcmp(argv[4], "default"))
        {
            PaToH_Initialize_Parameters(&args, PATOH_CONPART, PATOH_SUGPARAM_DEFAULT);
        }
        else
        {
            printf("An error occured during PaToH initilization. Please check your input parameters. \n");
            exit(1);
        }
    }
    else if (!strcmp(argv[2], "cutpart"))
    {
        if (!strcmp(argv[4], "quality"))
        {
            PaToH_Initialize_Parameters(&args, PATOH_CUTPART, PATOH_SUGPARAM_QUALITY);
        }
        else if (!strcmp(argv[4], "speed"))
        {
            PaToH_Initialize_Parameters(&args, PATOH_CUTPART, PATOH_SUGPARAM_SPEED);
        }
        else if (!strcmp(argv[4], "default"))
        {
            PaToH_Initialize_Parameters(&args, PATOH_CUTPART, PATOH_SUGPARAM_DEFAULT);
        }
        else
        {
            printf("An error occured during PaToH initilization. Please check your input parameters. \n");
            exit(1);
        }
    }
    else
    {
        printf("Please use either 'conpart' or 'cutpart' as the partition type.\n");
        exit(1);
    }

    args._k = atoi(argv[3]);
    args.seed = atoi(argv[5]);
    save_to_file = (atoi(argv[6]) == 0) ? false : true;

    PaToH_Alloc(&args, _c, _n, _nconst, cwghts, nwghts, xpins, pins);
    clock_t partitioning_begins = clock();
    PaToH_Part(&args, _c, _n, _nconst, 0, cwghts, nwghts, xpins, pins, NULL, partvec, partweights, &cut);
    clock_t partitioning_ends = clock();

    /***** Write the output files *****/

    if(save_to_file) {
      FILE *partvec_fp;
      partvec_fp = fopen(partvecFileName, "w");
      for (int i = 0; i < _c; i++)
      {
          fprintf(partvec_fp, "%d\n", partvec[i]);
      }
      fclose(partvec_fp);
    }

    PrintInfo(args._k, partweights, cut, _nconst, basename(argv[1]), partinfoFileName, nz, M, N);

    clock_t program_ends = clock();

    elapsedTime += (double)(program_ends - program_begins) / CLOCKS_PER_SEC;
    elapsed_time_patoh_part_only += (double)(partitioning_ends - partitioning_begins) / CLOCKS_PER_SEC;

    FILE *elapsedTime_fp;
    elapsedTime_fp = fopen(timeinfoFileName, "w");
    fprintf(elapsedTime_fp, "File Name: %s\n", basename(argv[1]));
    fprintf(elapsedTime_fp, "Non-Zero Count: %d\n", nz);
    fprintf(elapsedTime_fp, "Row and Column Count (M x N): %d and %d\n", M, N);
    fprintf(elapsedTime_fp, "Program took %.2f seconds to complete including IO operations.\n", elapsedTime);
    fprintf(elapsedTime_fp, "PaToH took %.2f seconds to partition the matrix.\n", elapsed_time_patoh_part_only);

    fclose(elapsedTime_fp);

    printf("Dataset: %s PartitionTime: %.2f\n", basename(argv[1]), elapsed_time_patoh_part_only);

    /***** Exit the program *****/
    free(cwghts);
    free(nwghts);
    free(partweights);
    free(partvec);
    PaToH_Free();
    exit(0);
}
