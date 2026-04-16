# hash_mult.h 函数分析文档

## 文件概述

`hash_mult.h` 是一个实现**稀疏矩阵乘法（SpGEMM - Sparse General Matrix Multiplication）**的头文件，使用**哈希表方法**来高效计算 `C = A × B`，其中 A、B、C 都是稀疏矩阵（CSR格式）。

### 核心特点
- 使用哈希表避免重复计算和排序
- 支持 OpenMP 多线程并行
- 支持 AVX2 向量化优化（SIMD）
- 两阶段执行：符号阶段（Symbolic Phase）和数值阶段（Numeric Phase）
- 支持多种数据类型：`int` 和 `long long int`

---

## 一、符号阶段（Symbolic Phase）函数

符号阶段的目标是**确定输出矩阵 C 每行的非零元素数量（NNZ）**，为后续内存分配做准备。

### 1. `symbolic_maxnnz_kernel<IT, NT>`
```cpp
template <class IT, class NT>
inline void symbolic_maxnnz_kernel(const IT *arpt, const IT *acol, 
                                    const IT *brpt, const IT *bcol, 
                                    BIN<IT, NT> &bin)
```
**功能**：计算每行的最大可能非零元素数量（上界估计）
- 通过简单累加 B 矩阵对应行的非零元素数量来估算
- 用于元数据输出（需要 `PRINT_METADATA_TO_FILE` 宏定义）

### 2. `hash_symbolic_kernel<IT, NT>`
```cpp
template <class IT, class NT>
inline void hash_symbolic_kernel(const IT *arpt, const IT *acol, 
                                  const IT *brpt, const IT *bcol,
                                  BIN<IT, NT> &bin, 
                                  const string& out_nnz_freq_outfile, 
                                  bool print_to_file = false)
```
**功能**：使用哈希表计算每行的精确非零元素数量
- **核心算法**：
  1. 为每行分配一个哈希表（大小根据 `bin_id` 动态调整）
  2. 遍历 A 的第 i 行，对每个非零元素 `A[i][j]`，遍历 B 的第 j 行
  3. 将 B 的列索引插入哈希表（使用线性探测解决冲突）
  4. 统计哈希表中的唯一列索引数量，即 `C[i]` 的非零元素数量
- **哈希函数**：`hash = (key * HASH_SCAL) & (ht_size - 1)`（位运算优化）
- 支持将结果输出到文件

### 3. `hash_symbolic_vec_kernel<NT>` (int 版本)
```cpp
template <class NT>
inline void hash_symbolic_vec_kernel(const int *arpt, const int *acol, 
                                      const int *brpt, const int *bcol, 
                                      BIN<int, NT> &bin)
```
**功能**：向量化版本的符号阶段（32位整数）
- 使用 **AVX2 SIMD 指令**（`__m256i`）同时处理 8 个元素（`VEC_LENGTH = 8`）
- 通过向量化比较和掩码操作加速哈希表查找
- 性能比标量版本快，但只支持 `int` 类型

### 4. `hash_symbolic_vec_kernel<NT>` (long long int 版本)
```cpp
template <class NT>
inline void hash_symbolic_vec_kernel(const long long int *arpt, ...)
```
**功能**：向量化版本的符号阶段（64位整数）
- 使用 AVX2 同时处理 4 个元素（`VEC_LENGTH_LONG = 4`）
- 因为 64 位整数在 256 位寄存器中只能放 4 个

### 5. `hash_symbolic<vectorProbing, IT, NT>`
```cpp
template <bool vectorProbing, class IT, class NT>
inline void hash_symbolic(const IT *arpt, const IT *acol, 
                          const IT *brpt, const IT *bcol, 
                          IT *crpt, BIN<IT, NT> &bin, 
                          const IT nrow, IT *nnz, 
                          string out_nnz_freq_outfile, 
                          bool first_pass = false)
```
**功能**：符号阶段的主入口函数
- 调用 `hash_symbolic_kernel` 或 `hash_symbolic_vec_kernel`（根据 `vectorProbing` 参数）
- 使用前缀和（scan）计算行指针 `crpt`
- 返回总非零元素数量 `nnz`

### 6. `hash_symbolic_topK<IT, NT>`
```cpp
template <class IT, class NT>
inline void hash_symbolic_topK(const IT *arpt, const IT *acol, ...)
```
**功能**：TopK 版本的符号阶段
- 用于只保留每行前 K 个最大值的场景
- 与标准版本类似，但为 TopK 优化

### 7. `hash_symbolic_kernel_with_inplace<IT, NT>`
```cpp
template <class IT, class NT>
inline void hash_symbolic_kernel_with_inplace(...)
```
**功能**：带原地更新优化的符号阶段
- 对于小行（非零元素数量 ≤ `inplace_cutoff`），直接预分配列索引数组
- 避免哈希表开销，提高小行的处理效率

---

## 二、数值阶段（Numeric Phase）函数

数值阶段的目标是**计算输出矩阵 C 的实际数值**，执行 `C[i][k] = Σ(A[i][j] * B[j][k])`。

### 1. `hash_numeric<sortOutput, IT, NT, MultiplyOperation, AddOperation>`
```cpp
template <bool sortOutput, typename IT, typename NT, 
          typename MultiplyOperation, typename AddOperation>
inline void hash_numeric(const IT *arpt, const IT *acol, const NT *aval,
                         const IT *brpt, const IT *bcol, const NT *bval,
                         const IT *crpt, IT *ccol, NT *cval, 
                         const BIN<IT, NT> &bin,
                         const MultiplyOperation multop, 
                         const AddOperation addop, IT cnnz)
```
**功能**：使用哈希表计算矩阵乘法的数值结果
- **核心算法**：
  1. 为每行分配哈希表（键：列索引，值：累加和）
  2. 遍历 A 的第 i 行，对每个 `A[i][j]`，遍历 B 的第 j 行
  3. 对每个 `B[j][k]`，计算 `A[i][j] * B[j][k]` 并累加到哈希表的 `C[i][k]` 位置
  4. 将哈希表内容写入输出矩阵 `ccol` 和 `cval`
- **参数说明**：
  - `sortOutput`：是否对输出列索引排序
  - `multop`：乘法操作（可以是普通乘法、自定义乘法等）
  - `addop`：加法操作（可以是普通加法、自定义累加等）

### 2. `hash_numeric_vec<sortOutput, NT, ...>` (int 版本)
```cpp
template <bool sortOutput, typename NT, ...>
inline void hash_numeric_vec(const int *arpt, ...)
```
**功能**：向量化版本的数值阶段（32位整数）
- 使用 AVX2 指令加速哈希表查找和更新
- 同时处理 8 个元素

### 3. `hash_numeric_vec<sortOutput, NT, ...>` (long long int 版本)
```cpp
template <bool sortOutput, typename NT, ...>
inline void hash_numeric_vec(const long long int *arpt, ...)
```
**功能**：向量化版本的数值阶段（64位整数）
- 使用 AVX2 同时处理 4 个元素

### 4. `hash_numeric_topk<sortOutput, IT, NT, ...>`
```cpp
template <bool sortOutput, typename IT, typename NT, ...>
inline void hash_numeric_topk(...)
```
**功能**：TopK 版本的数值阶段
- 只保留每行前 K 个最大值
- 使用 Jaccard 相似度或其他指标排序

### 5. `hash_numeric_with_inplace<sortOutput, IT, NT, ...>`
```cpp
template <bool sortOutput, typename IT, typename NT, ...>
inline void hash_numeric_with_inplace(...)
```
**功能**：带原地更新优化的数值阶段
- 小行（NNZ ≤ `inplace_cutoff`）：直接使用已排序的列索引数组，通过二分查找更新
- 大行：使用哈希表
- 减少内存分配和哈希表开销

### 6. `hash_numeric_with_inplace_V1<sortOutput, IT, NT, ...>`
```cpp
template <bool sortOutput, typename IT, typename NT, ...>
inline void hash_numeric_with_inplace_V1(...)
```
**功能**：原地更新 V1 版本
- 在符号阶段预填充列索引，数值阶段直接更新值
- 使用标志数组避免重复初始化

---

## 三、辅助函数

### 1. `sort_and_store_table2mat<sortOutput, IT, NT>`
```cpp
template <bool sortOutput, typename IT, typename NT>
inline void sort_and_store_table2mat(IT *ht_check, NT *ht_value, 
                                      IT *colids, NT *values, 
                                      IT nz, IT ht_size, IT offset, 
                                      IT row_id = -1)
```
**功能**：将哈希表内容排序并存储到输出矩阵
- 如果 `sortOutput = true`：对列索引排序后存储
- 如果 `sortOutput = false`：按哈希表顺序存储

### 2. `sort_and_store_table2mat_topK<sortOutput, IT, NT>`
**功能**：TopK 版本的排序和存储
- 按值排序，保留前 K 个

### 3. `sort_and_store_table2mat_topK_jaccard<sortOutput, IT, NT>`
**功能**：使用 Jaccard 相似度排序的 TopK 存储
- Jaccard 相似度 = `交集大小 / 并集大小`

### 4. `sort_and_store_table2mat_topK_cf<sortOutput, IT, NT>`
**功能**：使用 CF（可能是某种相似度指标）排序的 TopK 存储

### 5. `sanity_check<IT>`
```cpp
template <typename IT>
inline void sanity_check(const IT *crpt, IT cnnz, IT num_rows)
```
**功能**：检查输出矩阵的合法性
- 验证非零元素数量 > 0
- 验证行指针单调递增
- 验证行指针不超出范围

### 6. `sort_less<IT, NT>`
**功能**：用于 `std::sort` 的比较函数，按列索引升序排序

---

## 四、主要接口函数

### 1. `RowSpGEMM<vectorProbing, sortOutput, IT, NT, ...>`
```cpp
template <bool vectorProbing, bool sortOutput, 
          typename IT, typename NT, 
          typename MultiplyOperation, typename AddOperation>
void RowSpGEMM(const CSR<IT, NT> &a, const CSR<IT, NT> &b, CSR<IT, NT> &c,
               MultiplyOperation multop, AddOperation addop, 
               string out_nnz_freq_outfile, bool first_pass = false)
```
**功能**：**主要的稀疏矩阵乘法接口**
- **执行流程**：
  1. 初始化 BIN 结构（用于管理哈希表）
  2. 设置每行的最大 bin_id（根据预估的非零元素数量）
  3. 创建线程本地哈希表
  4. **符号阶段**：计算每行的非零元素数量
  5. 分配输出矩阵内存（`ccol` 和 `cval`）
  6. **数值阶段**：计算实际数值
- **参数**：
  - `vectorProbing`：是否使用向量化探测
  - `sortOutput`：是否对输出排序
  - `multop`/`addop`：自定义乘法和加法操作

### 2. `HashSpGEMMWithInplaceUpdate<vectorProbing, sortOutput, IT, NT, ...>`
**功能**：带原地更新优化的稀疏矩阵乘法
- 对小行使用原地更新，对大行使用哈希表

### 3. `HashSpGEMMWithInplaceUpdateV1<...>`
**功能**：原地更新 V1 版本
- 符号阶段预填充列索引

### 4. `HashSpGEMMtoCheckMemoryRequirement<...>`
**功能**：检查内存需求，不执行实际计算
- 用于评估内存使用量

### 5. `HashSpGEMMTopK<sortOutput, IT, NT, ...>`
**功能**：TopK 版本的稀疏矩阵乘法
- 只保留每行前 K 个最大值

---

## 五、技术细节

### 哈希表设计
- **哈希函数**：`hash = (key * HASH_SCAL) & (ht_size - 1)`
  - `HASH_SCAL = 107`（质数，减少冲突）
  - 使用位运算 `& (ht_size - 1)` 代替取模（要求 `ht_size` 是 2 的幂）
- **冲突解决**：线性探测（Linear Probing）
- **动态大小**：根据 `bin_id` 调整哈希表大小 `ht_size = MIN_HT_S << (bid - 1)`

### 并行化策略
- 使用 **OpenMP** 进行多线程并行
- 每个线程处理不同的行范围（`start_row` 到 `end_row`）
- 每个线程有独立的哈希表（避免锁竞争）

### 向量化优化
- 使用 **AVX2 SIMD 指令**（`__m256i`）
- 同时处理多个元素，减少循环次数
- 使用掩码操作（`_mm256_movemask_epi8`）快速查找匹配元素

### 内存优化
- **BIN 结构**：根据每行的非零元素数量动态分配哈希表大小
- **原地更新**：对小行直接使用已排序数组，避免哈希表开销
- **线程本地存储**：每个线程维护独立的哈希表，避免同步开销

---

## 六、使用示例

```cpp
// 基本用法
CSR<int, double> A, B, C;
// ... 初始化 A 和 B ...
RowSpGEMM<false, true>(A, B, C, 
                       [](double a, double b) { return a * b; },  // 乘法
                       [](double a, double b) { return a + b; },  // 加法
                       "output.txt", false);

// TopK 版本
HashSpGEMMTopK<true>(A, B, C, 
                     [](double a, double b) { return a * b; },
                     [](double a, double b) { return a + b; },
                     10,  // 保留前 10 个
                     "output.txt", false);
```

---

## 总结

`hash_mult.h` 实现了一个**高性能的稀疏矩阵乘法库**，主要特点：

1. **两阶段执行**：符号阶段确定结构，数值阶段计算值
2. **哈希表优化**：避免重复计算和排序开销
3. **多线程并行**：OpenMP 并行处理
4. **向量化加速**：AVX2 SIMD 指令优化
5. **内存优化**：动态哈希表大小、原地更新等
6. **灵活接口**：支持自定义操作、TopK、排序等

适用于大规模稀疏矩阵乘法计算，如图分析、机器学习、科学计算等领域。

