/* ************************************************************************
 * Copyright (c) 2026 Advanced Micro Devices, Inc.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * ************************************************************************ */

#include "aoclsparse.h"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

struct Entry
{
    aoclsparse_int row;
    aoclsparse_int col;
    double         val;
};

static void check_error(aoclsparse_status status, const char *api_name)
{
    if(status != aoclsparse_status_success)
    {
        std::cerr << "ERROR in " << api_name << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

static std::string matrix_name_from_path(const std::string &path)
{
    std::size_t pos = path.find_last_of("/\\");
    std::string name = (pos == std::string::npos) ? path : path.substr(pos + 1);
    std::size_t dot = name.find_last_of('.');
    if(dot != std::string::npos)
        name = name.substr(0, dot);
    return name;
}

static bool read_matrix_market_to_csr(const std::string               &filename,
                                      aoclsparse_int                  &m,
                                      aoclsparse_int                  &n,
                                      std::vector<aoclsparse_int>     &row_ptr,
                                      std::vector<aoclsparse_int>     &col_ind,
                                      std::vector<double>             &vals)
{
    std::ifstream fin(filename);
    if(!fin.is_open())
        return false;

    std::string line;
    if(!std::getline(fin, line))
        return false;

    bool symmetric = (line.find("symmetric") != std::string::npos)
                     || (line.find("hermitian") != std::string::npos);
    bool pattern = (line.find("pattern") != std::string::npos);

    do
    {
        if(!std::getline(fin, line))
            return false;
    } while(!line.empty() && line[0] == '%');

    std::istringstream size_ss(line);
    aoclsparse_int      nnz_file = 0;
    size_ss >> m >> n >> nnz_file;
    if(size_ss.fail() || m <= 0 || n <= 0 || nnz_file < 0)
        return false;

    std::vector<Entry> entries;
    entries.reserve(static_cast<std::size_t>(symmetric ? 2 * nnz_file : nnz_file));

    for(aoclsparse_int i = 0; i < nnz_file; ++i)
    {
        if(!std::getline(fin, line))
            return false;
        if(line.empty() || line[0] == '%')
        {
            --i;
            continue;
        }

        std::istringstream iss(line);
        aoclsparse_int r1 = 0, c1 = 0;
        double         v = 1.0;
        iss >> r1 >> c1;
        if(iss.fail())
            return false;
        if(!pattern)
            iss >> v;

        aoclsparse_int r = r1 - 1;
        aoclsparse_int c = c1 - 1;
        entries.push_back({r, c, v});
        if(symmetric && r != c)
            entries.push_back({c, r, v});
    }

    std::sort(entries.begin(),
              entries.end(),
              [](const Entry &a, const Entry &b) {
                  if(a.row != b.row)
                      return a.row < b.row;
                  return a.col < b.col;
              });

    // Merge duplicate coordinates by summation.
    std::vector<Entry> merged;
    merged.reserve(entries.size());
    for(const Entry &e : entries)
    {
        if(!merged.empty() && merged.back().row == e.row && merged.back().col == e.col)
            merged.back().val += e.val;
        else
            merged.push_back(e);
    }

    row_ptr.assign(static_cast<std::size_t>(m + 1), 0);
    for(const Entry &e : merged)
        row_ptr[static_cast<std::size_t>(e.row + 1)]++;
    for(aoclsparse_int i = 0; i < m; ++i)
        row_ptr[static_cast<std::size_t>(i + 1)] += row_ptr[static_cast<std::size_t>(i)];

    col_ind.resize(merged.size());
    vals.resize(merged.size());
    for(std::size_t i = 0; i < merged.size(); ++i)
    {
        col_ind[i] = merged[i].col;
        vals[i] = merged[i].val;
    }

    return true;
}

struct SpgemmTwoStageTime
{
    double nnz_count_ms;
    double finalize_ms;
    double total_ms;
};

static SpgemmTwoStageTime run_spgemm_two_stage_once(aoclsparse_mat_descr descrA,
                                                    aoclsparse_matrix    csrA,
                                                    aoclsparse_mat_descr descrB,
                                                    aoclsparse_matrix    csrB)
{
    const aoclsparse_operation transA = aoclsparse_operation_none;
    const aoclsparse_operation transB = aoclsparse_operation_none;
    aoclsparse_matrix          csrC = nullptr;

    auto t0 = std::chrono::high_resolution_clock::now();
    aoclsparse_status status = aoclsparse_dcsr2m(transA,
                                                 descrA,
                                                 csrA,
                                                 transB,
                                                 descrB,
                                                 csrB,
                                                 aoclsparse_stage_nnz_count,
                                                 &csrC);
    auto t1 = std::chrono::high_resolution_clock::now();
    check_error(status, "aoclsparse_scsr2m(stage_nnz_count)");

    status = aoclsparse_dcsr2m(
        transA, descrA, csrA, transB, descrB, csrB, aoclsparse_stage_finalize, &csrC);
    auto t2 = std::chrono::high_resolution_clock::now();
    check_error(status, "aoclsparse_scsr2m(stage_finalize)");

    aoclsparse_destroy(&csrC);

    std::chrono::duration<double, std::milli> nnz_ms = t1 - t0;
    std::chrono::duration<double, std::milli> fin_ms = t2 - t1;
    return {nnz_ms.count(), fin_ms.count(), nnz_ms.count() + fin_ms.count()};
}

int main(int argc, char **argv)
{
    if(argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <matrix.mtx>" << std::endl;
        return EXIT_FAILURE;
    }

    const std::string mtx_path(argv[1]);
    const std::string mtx_name = matrix_name_from_path(mtx_path);

    aoclsparse_int              m = 0, n = 0;
    std::vector<aoclsparse_int> row_ptr, col_ind;
    std::vector<double>         vals;
    if(!read_matrix_market_to_csr(mtx_path, m, n, row_ptr, col_ind, vals))
    {
        std::cerr << "Failed to read matrix file: " << mtx_path << std::endl;
        return EXIT_FAILURE;
    }

    aoclsparse_matrix    csrA = nullptr;
    aoclsparse_matrix    csrB = nullptr;
    aoclsparse_mat_descr descrA;
    aoclsparse_mat_descr descrB;
    aoclsparse_index_base base = aoclsparse_index_base_zero;

    check_error(aoclsparse_create_mat_descr(&descrA), "aoclsparse_create_mat_descr(A)");
    check_error(aoclsparse_create_mat_descr(&descrB), "aoclsparse_create_mat_descr(B)");

    // Benchmark A*A (square matrix assumed by the benchmark usage).
    check_error(aoclsparse_create_dcsr(
                    &csrA,
                    base,
                    m,
                    n,
                    static_cast<aoclsparse_int>(vals.size()),
                    row_ptr.data(),
                    col_ind.data(),
                    vals.data()),
                "aoclsparse_create_dcsr(A)");
    check_error(aoclsparse_create_dcsr(
                    &csrB,
                    base,
                    m,
                    n,
                    static_cast<aoclsparse_int>(vals.size()),
                    row_ptr.data(),
                    col_ind.data(),
                    vals.data()),
                "aoclsparse_create_dcsr(B)");

    // Warmup: run both two-stage steps once.
    (void)run_spgemm_two_stage_once(descrA, csrA, descrB, csrB);

    // Timed runs
    constexpr int iters = 10;
    double        nnz_count_ms = 0.0;
    double        finalize_ms = 0.0;
    double        total_ms = 0.0;
    for(int i = 0; i < iters; ++i)
    {
        SpgemmTwoStageTime t = run_spgemm_two_stage_once(descrA, csrA, descrB, csrB);
        nnz_count_ms += t.nnz_count_ms;
        finalize_ms += t.finalize_ms;
        total_ms += t.total_ms;
    }

    const double avg_nnz_count_ms = nnz_count_ms / static_cast<double>(iters);
    const double avg_finalize_ms = finalize_ms / static_cast<double>(iters);
    const double avg_total_ms = total_ms / static_cast<double>(iters);
    std::cout << mtx_name << ", " << std::fixed << std::setprecision(6) << avg_nnz_count_ms << ", "
              << avg_finalize_ms << ", " << avg_total_ms << std::endl;

    aoclsparse_destroy_mat_descr(descrA);
    aoclsparse_destroy_mat_descr(descrB);
    aoclsparse_destroy(&csrA);
    aoclsparse_destroy(&csrB);

    return EXIT_SUCCESS;
}
