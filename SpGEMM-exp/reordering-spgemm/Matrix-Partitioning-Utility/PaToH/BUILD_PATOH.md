# PaToH 从源代码编译指南

## 问题分析

`kmer_V1r` 是一个非常大的矩阵（214M x 214M，4.2GB），使用 `quality` 模式进行划分需要很长时间：
- k64 版本：约 19 小时（68632 秒）
- k56 版本：可能更长，当前已运行 24+ 小时

## 方案1：从源代码编译 PaToH（不推荐，因为没有源代码提供）

### 步骤1：下载 PaToH 

PaToH 的官方下载地址：
```
https://faculty.cc.gatech.edu/~umit/software.html
```

### 步骤4：复制到项目目录

```bash
# 复制库文件
cp build/Linux-x86_64/libpatoh.a /data/lsl/SpGEMM/SpGEMM-exp/reordering-spgemm/Matrix-Partitioning-Utility/PaToH/libpatoh_linux.a

# 复制头文件（如果需要）
cp build/Linux-x86_64/patoh.h /data/lsl/SpGEMM/SpGEMM-exp/reordering-spgemm/Matrix-Partitioning-Utility/PaToH/patoh.h
```

### 步骤5：重新编译主程序

```bash
cd /data/lsl/SpGEMM/SpGEMM-exp/reordering-spgemm/Matrix-Partitioning-Utility/PaToH
gcc main.c libpatoh_linux.a -lm -O3 -o a.out
```

## 方案2：优化当前运行（临时方案）

### 选项A：使用 speed 或 default 模式

对于超大规模矩阵，`quality` 模式可能过于耗时。可以尝试：

```bash
# 使用 default 模式（平衡质量和速度）
./a.out "${DATA_PATH}/kmer_V1r/kmer_V1r.mtx" cutpart 56 default 1 1

# 或使用 speed 模式（最快，但质量可能略低）
./a.out "${DATA_PATH}/kmer_V1r/kmer_V1r.mtx" cutpart 56 speed 1 1
```

### 选项B：检查当前进程状态

当前进程仍在运行，可以监控其进度：

```bash
# 查看进程状态
ps aux | grep kmer_V1r

# 查看输出文件是否在更新
ls -lth output/PaToH_kmer_V1r.mtx_cutpart_k56* 2>/dev/null

# 如果文件存在，查看最后几行
tail -f output/PaToH_kmer_V1r.mtx_cutpart_k56_quality_s1_timeinfo.txt
```

### 选项C：使用 k64 的结果（如果可接受）

如果 k64 的结果已经足够好，可以考虑：
1. 直接使用 k64 的划分结果
2. 或者将 k64 的结果转换为 k56（需要额外的后处理）

## 方案3：使用 METIS 作为替代

如果 PaToH 太慢，可以考虑使用 METIS 进行划分：

```bash
cd /data/lsl/SpGEMM/SpGEMM-exp/reordering-spgemm/Matrix-Partitioning-Utility/METIS
./run -i "${DATA_PATH}/kmer_V1r/kmer_V1r.mtx" -k 56 -o edge-cut -t -s
```

METIS 通常比 PaToH 更快，但划分质量可能略有不同。

## 性能优化建议

1. **使用多线程版本**（如果 PaToH 支持）
2. **增加系统内存**（当前使用约 85GB）
3. **使用更快的存储**（SSD vs HDD）
4. **调整 PaToH 内部参数**（需要修改源代码）

## 注意事项

- 从源代码编译可能需要一些依赖库
- 编译选项（如 `-O3`）可能影响性能和正确性
- 建议在测试数据集上先验证编译后的库是否正常工作

## 参考链接

- PaToH 官方网站：https://faculty.cc.gatech.edu/~umit/software.html
- PaToH 文档：查看下载包中的 README 文件

