#ifndef _IO_SPGEMM_H
#define _IO_SPGEMM_H

#include <cstdio>
#include <cassert>
#include <type_traits>
#include "Triple.h"
#include "CSC.h"

#define READBUFFER (512 * 1024 * 1024)  // in MB

template <typename IT, typename NT>
int ReadBinary(string filename, CSC<IT,NT> & csc)
{
    FILE * f = fopen(filename.c_str(), "r");
    if(!f)
    {
        cerr << "Problem reading binary input file" << filename << endl;
        return -1;
    }
    IT m,n,nnz;
    fread(&m, sizeof(IT), 1, f);
    fread(&n, sizeof(IT), 1, f);
    fread(&nnz, sizeof(IT), 1, f);
    
    if (m <= 0 || n <= 0 || nnz <= 0)
    {
        cerr << "Problem with matrix size in binary input file" << filename << endl;
        return -1;
    }
    double start = omp_get_wtime( );
    cout << "Reading matrix with dimensions: "<< m << "-by-" << n <<" having "<< nnz << " nonzeros" << endl;
    
    IT * rowindices = new IT[nnz];
    IT * colindices = new IT[nnz];
    NT * vals = new NT[nnz];
    
    size_t rows = fread(rowindices, sizeof(IT), nnz, f);
    size_t cols = fread(colindices, sizeof(IT), nnz, f);
    size_t nums = fread(vals, sizeof(NT), nnz, f);
    
    if(rows != nnz || cols != nnz || nums != nnz)
    {
        cerr << "Problem with FREAD, aborting... " << endl;
        return -1;
    }
    double end = omp_get_wtime( );
    // printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);
    printf("Converting matrix data from binary to COO fromat takes %.16g seconds.\n", end - start);

    fclose(f);
    
    csc = *(new CSC<IT,NT>(rowindices, colindices, vals , nnz, m, n));
    
    delete [] rowindices;
    delete [] colindices;
    delete [] vals;
    return 1;
}

template <typename IT, typename NT>
int ReadASCII(string filename, CSC<IT,NT> & csc, IT cluster_sz = 0)
{
    bool isSymmetric = false;
    double start = omp_get_wtime( );
    
    FILE *fid = fopen(filename.c_str(), "r");
    if (!fid) {
        cerr << "Problem reading ASCII input file: " << filename << endl;
        return -1;
    }
    
    // Read header lines (starting with %)
    char line[256];
    bool found_size_line = false;
    while (fgets(line, sizeof(line), fid) != NULL) {
        // Skip blank lines
        int len = strlen(line);
        bool is_blank = true;
        for (int i = 0; i < len; i++) {
            if (line[i] != ' ' && line[i] != '\t' && line[i] != '\n' && line[i] != '\r') {
                is_blank = false;
                break;
            }
        }
        if (is_blank) continue;
        
        if (line[0] == '%') {
            if (strstr(line, "symmetric")) {
                isSymmetric = true;
                cout << "Matrix is symmetric" << endl;
            }
        } else {
            // This should be the size line
            found_size_line = true;
            break;
        }
    }
    
    if (!found_size_line) {
        cerr << "Error: Could not find size line in " << filename << endl;
        fclose(fid);
        return -1;
    }
    
    // Parse size line
    IT m, n, nnz;
    int items_read = 0;
    
    if constexpr(std::is_same<IT, int>::value) {
        items_read = sscanf(line, "%d %d %d", &m, &n, &nnz);
    } else if constexpr(std::is_same<IT, long>::value) {
        items_read = sscanf(line, "%ld %ld %ld", &m, &n, &nnz);
    } else {
        // For long long or int64_t, use %lld
        items_read = sscanf(line, "%lld %lld %lld", &m, &n, &nnz);
    }
    
    if (items_read != 3) {
        cerr << "Error: Failed to parse size line from " << filename << endl;
        cerr << "Line content: [" << line << "]" << endl;
        cerr << "Expected 3 values (m n nnz), got " << items_read << endl;
        cerr << "IT type is: ";
        if constexpr(std::is_same<IT, int>::value) {
            cerr << "int" << endl;
        } else if constexpr(std::is_same<IT, long>::value) {
            cerr << "long (using %ld format)" << endl;
        } else {
            cerr << "long long / int64_t (using %lld format)" << endl;
        }
        fclose(fid);
        return -1;
    }
    
    IT original_nnz = nnz;
    IT max_nnz = nnz;
    if (isSymmetric) {
        max_nnz = nnz * 2;  // maximum possible nnz (if no diagonal elements)
    }
    
    Triple<IT,NT> * triples = new Triple<IT,NT>[max_nnz];
    IT cnz = 0;  // current number of nonzeros
    
    // Read triplets using fscanf
    for (IT i = 0; i < original_nnz; i++) {
        IT row_id, col_id;
        double V = 1.0;  // default value for binary mode (row, col only)
        int items_read;
        
        // Read row and column first
        if constexpr(std::is_same<IT, int>::value) {
            items_read = fscanf(fid, "%d %d", &row_id, &col_id);
        } else if constexpr(std::is_same<IT, long>::value) {
            items_read = fscanf(fid, "%ld %ld", &row_id, &col_id);
        } else {
            // For long long or int64_t, use %lld
            items_read = fscanf(fid, "%lld %lld", &row_id, &col_id);
        }
        
        assert(items_read == 2 && "Failed to read row and column");
        
        // Try to read value, if not available (binary mode), use default 1.0
        // Check if there's more data on the line (not just whitespace/newline)
        char next_char = fgetc(fid);
        if (next_char != '\n' && next_char != EOF && next_char != '\r') {
            ungetc(next_char, fid);
            // Try to read a double value
            if (fscanf(fid, "%lf", &V) != 1) {
                V = 1.0;  // default value if read fails
            }
            // Skip remaining whitespace to end of line
            while ((next_char = fgetc(fid)) != '\n' && next_char != EOF && next_char != '\r');
        }
        
        IT row_0based = (IT)row_id - 1;  // convert to 0-based
        IT col_0based = (IT)col_id - 1;  // convert to 0-based
        
        triples[cnz].row = row_0based;
        triples[cnz].col = col_0based;
        triples[cnz].val = (NT)V;
        
        if (isSymmetric) {
            if (row_0based != col_0based) {
                // Non-diagonal element: add symmetric entry
                cnz++;
                triples[cnz].col = row_0based;
                triples[cnz].row = col_0based;
                triples[cnz].val = (NT)V;
            }
            // Diagonal element: don't duplicate, just use the single entry
        }
        ++cnz;
    }
    
    // Update nnz to actual number of nonzeros read (including symmetric entries)
    nnz = cnz;
    fclose(fid);
    
    double end = omp_get_wtime( );
    printf("Converting matrix data from ASCII to COO format: %.16g seconds\n", end - start);
    printf("Input Matrix: Rows = %d, Columns= %d, nnz = %d, cluster_sz = %d\n", m, n, nnz, cluster_sz);

    cout << "Converting to csc ... " << endl << endl;
    if(cluster_sz != 0) {
      if(m%cluster_sz != 0) {
        m = ((m / cluster_sz) + 1) * cluster_sz;
        n = m;
      }
    }
    csc= *(new CSC<IT,NT>(triples, nnz, m, n));
    csc.totalcols = n;
    delete [] triples;
    return 1;
}

#endif
