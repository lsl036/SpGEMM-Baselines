# Kokkos-Kernels SpGEMM 测试说明

## 0. 为什么 build 里没有生成 Makefile？

**原因**：kokkos-kernels 是 Kokkos 的上层库，CMake 配置时会执行 `find_package(Kokkos REQUIRED)`。**如果找不到 Kokkos，CMake 会直接报错退出，不会生成任何 Makefile**。

因此必须先 **先编译并安装 Kokkos**，再在配置 kokkos-kernels 时 **指定 Kokkos 的路径**（见下文）。

---

## 1. 正确编译方式（两步）

### 步骤 1：先编译并安装 Kokkos

若你已有 Kokkos 源码（例如在 `SpGEMM/kokkos/`）：

```bash
cd /path/to/kokkos
mkdir -p build && cd build
cmake .. \
  -DCMAKE_INSTALL_PREFIX=/path/to/kokkos/install \
  -DKokkos_ENABLE_SERIAL=On \
  -DKokkos_ENABLE_OPENMP=On
make -j
make install
```

`/path/to/kokkos/install` 可改成你想要的安装目录；若用 CUDA/HIP 等，需加上对应 `-DKokkos_ENABLE_*`。

### 步骤 2：再配置并编译 kokkos-kernels（必须指定 Kokkos）

**当前环境**：Kokkos 已安装在 **`/usr/local`**（构建目录在 `SpGEMM/kokkos/builddir`）。任选其一即可。

**若你刚才已经 cmake 成功但 Makefile 在源码目录**（`build/` 里没有 Makefile）：多半是命令行里多了空格被 CMake 当成路径，导致构建目录错了。你可以**先直接编译**：

```bash
cd /data/lsl/SpGEMM/kokkos-kernels
make -j
```

Makefile 在 `kokkos-kernels/` 根目录时，在根目录执行 `make -j` 即可；可执行文件会生成在 `kokkos-kernels/perf_test/sparse/sparse_spgemm` 等相对路径下。

**下次希望 Makefile 在 `build/` 里**，请用 **`-S` / `-B` 显式指定源码和构建目录**（避免续行空格被当路径）：

```bash
cd /data/lsl/SpGEMM/kokkos-kernels
rm -rf build && mkdir -p build
# 方式 A：用已安装的 Kokkos
cmake -S . -B build -DKokkos_ROOT=/usr/local -DKokkosKernels_INST_DOUBLE=On -DCMAKE_INSTALL_PREFIX=./ -DKokkosKernels_ENABLE_TESTS_AND_PERFSUITE=ON -DKokkosKernels_ENABLE_TESTS=ON -DKokkosKernels_ENABLE_PERFTESTS=ON
# 方式 B：用 Kokkos 构建目录
# cmake -S . -B build -DKokkos_DIR=/data/lsl/SpGEMM/kokkos/builddir -DKokkosKernels_ENABLE_PERFTESTS=ON ...
cd build && make -j
```

**要点**：

- **必须加 `-DKokkosKernels_ENABLE_PERFTESTS=ON`**，否则不会编译 `perf_test`（也就没有 `sparse_spgemm` 等性能测试可执行文件）。
- **`-DKokkos_ROOT=/usr/local`** 指向已安装的 Kokkos。若只用构建目录，用 **`-DKokkos_DIR=/data/lsl/SpGEMM/kokkos/builddir`**。
- 使用 **`-S . -B build`** 可以保证构建文件一定在 `build/`，不会因续行空格导致生成到源码目录。
- 若 CMake 生成的是 Ninja，则在 `build/` 里执行 **`ninja`** 代替 `make -j`。

---

## 2. SpGEMM 可执行文件位置

编译成功后，SpGEMM 相关的 **performance test** 可执行文件在：

| 可执行文件 | 路径（out-of-source 在 build/ 时） | 路径（in-source 时） |
|------------|--------------------------------------|----------------------|
| **sparse_spgemm** | `build/perf_test/sparse/sparse_spgemm` | `perf_test/sparse/sparse_spgemm` |
| **sparse_spgemm_jacobi** | `build/perf_test/sparse/sparse_spgemm_jacobi` | `perf_test/sparse/sparse_spgemm_jacobi` |

查找命令：

```bash
find . -name 'sparse_spgemm*' -type f
```

---

## 3. SpGEMM 测试/实现代码位置

### 3.1 Performance Test（跑性能、读矩阵、计时）

| 用途 | 路径 |
|------|------|
| SpGEMM 性能测试主程序 | **`perf_test/sparse/KokkosSparse_spgemm.cpp`** |
| SpGEMM Jacobi 性能测试 | **`perf_test/sparse/KokkosSparse_spgemm_jacobi.cpp`** |

- 入口：`run_spgemm()`（根据 execution space 模板实例化）
- 命令行：`--amtx`、`--bmtx` 指定矩阵（支持 `.mtx`），不写 `--bmtx` 则算 C=AxA
- 读矩阵：`KokkosSparse::Impl::read_kokkos_crst_matrix()`，在 **`sparse/src/KokkosSparse_IOUtils.hpp`**（支持 `.mtx` / `.mm`）

### 3.2 SpGEMM 内核与 API

| 用途 | 路径 |
|------|------|
| SpGEMM 接口与 handle | `sparse/src/KokkosSparse_spgemm.hpp`、`KokkosSparse_spgemm_handle.hpp` |
| 实现（symbolic/numeric 等） | `sparse/impl/KokkosSparse_spgemm_impl*.hpp`、`KokkosSparse_spgemm_*_spec.hpp` |
| Unit Test | `sparse/unit_test/Test_Sparse_spgemm.hpp` |
| 示例 | `example/wiki/sparse/KokkosSparse_wiki_spgemm.cpp` |

---

## 4. 运行 SpGEMM 性能测试

```bash
cd build
# C = A * A（只给 A）
./perf_test/sparse/sparse_spgemm --amtx /path/to/matrix.mtx --repeat 5

# C = A * B
./perf_test/sparse/sparse_spgemm --amtx /path/to/A.mtx --bmtx /path/to/B.mtx --repeat 5

# 选算法、写结果、检查正确性
./perf_test/sparse/sparse_spgemm --amtx A.mtx --algorithm KKDEFAULT --verbose --checkoutput --cmtx C.mtx
```

`--algorithm` 可选：`DEFAULT` / `KKDEFAULT` / `KKSPGEMM`、`KKMEM`、`KKDENSE` 等（见 `print_options()` 或 `--help` 输出）。

---

## 5. 性能测试很慢时检查（优化 / 并行）

若 `sparse_spgemm` 跑得很慢，按下面逐项确认：

| 检查项 | 说明 |
|--------|------|
| **1) 构建类型** | 必须用 **Release**。配置时加 **`-DCMAKE_BUILD_TYPE=Release`**（Kokkos 和 kokkos-kernels 都要）。Debug 或未设置会慢一个数量级以上。 |
| **2) 是否在用 OpenMP** | 运行时**必须**加 **`--openmp N`**（N=核心数）。不加会走 Serial，且若 Serial 未 enable 会报错。终端里应看到 `Running on OpenMP backend.` 和 `thread_pool_topology[ 1 x N x 1 ]`。 |
| **3) OpenMP 绑定** | 建议设置 **`export OMP_PROC_BIND=spread OMP_PLACES=threads`**，否则 Kokkos 会打 WARNING，并行效率可能变差。`run_spgemm_benchmark.sh` 已自动设置。 |
| **4) 算法** | 默认已是 **KKSPGEMM**。大矩阵可尝试 **`--algorithm KKDENSE`** 或 **KKMEM** 做对比。 |
| **5) 线程数** | `--openmp N` 的 N 建议设为实际可用 CPU 核心数（如 28）。 |

**推荐 cmake 配置（Release + 显式 OpenMP）：**
```bash
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DKokkos_ROOT=/usr/local \
  -DKokkosKernels_INST_DOUBLE=On \
  -DKokkosKernels_ENABLE_PERFTESTS=ON \
  ...
cd build && make -j
```

**推荐运行方式：**
```bash
export OMP_PROC_BIND=spread OMP_PLACES=threads
./perf_test/sparse/sparse_spgemm --amtx matrix.mtx --repeat 5 --openmp 28
```

---

## 6. 若 build 里没有 Makefile 或没有可执行文件

- **build 里根本没有 Makefile**：说明 CMake 配置阶段就失败了，多半是 **没找到 Kokkos**。请按 **第 0、1 节** 先安装 Kokkos，并用 `-DKokkos_ROOT=...` 或 `-DKokkos_DIR=...` 再跑一次 `cmake ..`，确认终端里有 “Build files have been written to ...” 再执行 `make` 或 `ninja`。
- **perf_test/sparse 没有生成 sparse_spgemm**：需要同时满足两点：
  1. **`KokkosKernels_ENABLE_PERFTESTS=ON`**  
     否则不会进入 `perf_test/`，整个 `perf_test/sparse/` 都不会被编译。  
  2. **Sparse 组件打开**  
     `perf_test/sparse/` 由 component “sparse” 控制。默认 **`KokkosKernels_ENABLE_ALL_COMPONENTS=ON`** 会打开包括 sparse 在内的所有组件；若你曾把 `ENABLE_ALL_COMPONENTS` 设为 `OFF`，需要显式加上 **`-DKokkosKernels_ENABLE_COMPONENT_SPARSE=ON`**。  

  **推荐完整配置（保证有 sparse_spgemm）：**
  ```bash
  cmake -S . -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DKokkos_ROOT=/usr/local \
    -DKokkosKernels_INST_DOUBLE=On \
    -DKokkosKernels_ENABLE_PERFTESTS=ON \
    -DKokkosKernels_ENABLE_TESTS_AND_PERFSUITE=ON \
    -DKokkosKernels_ENABLE_TESTS=ON
  cd build && make -j
  ```
  可执行文件在 **`build/perf_test/sparse/sparse_spgemm`**（out-of-source）或 **`perf_test/sparse/sparse_spgemm`**（in-source 时在源码根目录下）。用 `find . -name sparse_spgemm -type f` 可确认位置。
- **运行时报错 "Serial is not enabled"**：当前安装的 Kokkos 只开了 OpenMP、没开 Serial，程序在单线程或未请求并行时会回退到 Serial 导致报错。需要 **重新配置并安装 Kokkos**，加上 **`-DKokkos_ENABLE_SERIAL=On`**（与 OpenMP 同时开），再 `make install`，然后重新编译 kokkos-kernels。
- **链接报错 undefined reference to ... Kokkos::Serial ...（KokkosKernels_sparse_serial 等）**：单元测试需要 **Serial 后端的 ETI**（显式模板实例化），但 libkokkoskernels 里没有。两种办法二选一：  
  **办法 A（只跑 SpGEMM 性能测试）**：关闭单元测试，只开性能测试，不编 Serial 单元测试即可：  
  ```bash
  rm -rf build && mkdir build && cd build
  cmake -S .. -B . -DCMAKE_BUILD_TYPE=Release -DKokkos_ROOT=/usr/local \
    -DKokkosKernels_INST_DOUBLE=On -DKokkosKernels_ENABLE_PERFTESTS=ON \
    -DKokkosKernels_ENABLE_TESTS=OFF
  make -j
  ```  
  **办法 B（要跑单元测试）**：保证 Kokkos 已带 Serial，再**清空构建目录**后重新配置，让 kokkos-kernels 打开 Serial ETI 并重编库：  
  ```bash
  # 确认 Kokkos 有 Serial：grep Kokkos_ENABLE_SERIAL /usr/local/lib/cmake/Kokkos/KokkosConfig.cmake
  rm -rf build && mkdir build && cd build
  cmake -S .. -B . -DCMAKE_BUILD_TYPE=Release -DKokkos_ROOT=/usr/local \
    -DKokkosKernels_INST_DOUBLE=On -DKokkosKernels_ENABLE_PERFTESTS=ON \
    -DKokkosKernels_ENABLE_TESTS=ON -DKokkosKernels_INST_EXECSPACE_SERIAL=ON
  make -j
  ```
- **有 Makefile 但没有可执行文件（且已开 PERFTESTS）**：看完整编译错误：`make -j 2>&1 | tee build.log`；并确认 Kokkos 与 kokkos-kernels 的配置一致（如都开启 Serial/OpenMP 等）。
