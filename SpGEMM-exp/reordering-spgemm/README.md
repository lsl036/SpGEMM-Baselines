# Reordering SpGEMM

This artifact integrates a set of tools to perform reordering of sparse
matrices using the ten algorithms listed in the following table. It primarily
serves as a preprocessing step of our [Cluster-wise SpGEMM](https://github.com/PASSIONLab/clusterwise-spgemm) evaluation.

For more information, please read our preprint paper ["Improving SpGEMM Performance Through Matrix Reordering and Cluster-wise Computation"](https://arxiv.org/abs/2507.21253) on arXiv.

| **Algorithm**                    | **Description** | **Original Repo**                                                                                   |
|----------------------------------|-----------------|-----------------------------------------------------------------------------------------------------|
| Original                         | Original input order |  |
| Random                           | Random shuffle |  |
| Rev. Cuthill–McKee (RCM) [1]     | Bandwidth reduction via BFS | [libmtx](https://github.com/simulahpc/libmtx)                                                       |
| Aprx. minimum degree (AMD) [2]   | Greedy strategy to reduce fill | [SparseBase](https://github.com/sparcityeu/SparseBase) |
| Nested dissection (ND) [3]       | Recursive divide-and-conquer to reduce fill | [libmtx](https://github.com/simulahpc/libmtx)                                                       |
| Graph partitioning (GP) [4]      | METIS using edge-cut objective | [METIS](https://karypis.github.io/glaros/files/sw/metis/metis-5.1.0.tar.gz), [Matrix-Partitioning-Utility](https://github.com/sparcityeu/Matrix-Partitioning-Utility) |
| Hypergraph partitioning (HP) [5] | PaToH using cut-net metric | [Matrix-Partitioning-Utility](https://github.com/sparcityeu/Matrix-Partitioning-Utility)            |
| Gray code ordering [6]           | Splitting sparse and dense rows | [SparseBase](https://github.com/sparcityeu/SparseBase) |
| Rabbit order [7]                 | Hierarchical community-based reordering | [rabbit_order](https://github.com/araij/rabbit_order) |
| Degree order                     | Reorder in descending order of degrees | [SparseBase](https://github.com/sparcityeu/SparseBase) |
| Slash-burn method (SB) [8]       | Recursively split rows into hubs and spokes | [SparseBase](https://github.com/sparcityeu/SparseBase) |

## Requirements
### Hardware
This artifact has no specific hardware requirements and is expected to operate on any machine with at least 192 GB memory and the required software stack installed. However, for consistency, this artifact should be executed on the same platform intended for evaluating [clusterwise-spgemm](https://github.com/PaSSIONLab/clusterwise-spgemm).

### Software
As this artifact integrates a set of tools implementing ten reordering algorithms across different libraries, it requires several software dependencies, including:
1. **RCM & ND** require the installation of [libmtx-v0.5.0](https://github.com/simulahpc/libmtx), which further depends on standard GNU build tools, including Make, Autoconf, Automake, and Libtool. The source code must be compiled with a C99-compliant compiler (e.g., GCC or Clang) and requires BLAS and MPI libraries for core linear algebra and parallel communication functionality. Optional support for OpenMP enables multi-threaded execution.
2. **GP** requires a C++11-compliant compiler (e.g., GCC or Clang), standard Unix development tools (e.g., GNU Make), and the [METIS 5.1.0](https://karypis.github.io/glaros/files/sw/metis/metis-5.1.0.tar.gz) graph partitioning library.
3. **HP** requires a standard C compiler (e.g., GCC) and a precompiled PaToH static library (provided in the repository), with no additional external dependencies.
4. **AMD, Gray, Degree & SlashBurn** are generated using [SparseBase-v0.3.1](https://github.com/sparcityeu/SparseBase). SparseBase requires CMake (version 3.12 or later) and a C++17-compliant compiler (such as GCC or Clang). To enable AMD support, [SuiteSparse](https://github.com/DrTimothyAldenDavis/SuiteSparse) must be installed beforehand. Further installation details are available in the official [SparseBase documentation](https://sparcityeu.github.io/SparseBase/).
5. **Rabbit** is generated using the [official GitHub repository](https://github.com/araij/rabbit_order). It requires G++ (≥v4.9.2) with C++11 support and depends on the Boost C++ Libraries (≥v1.58.0) for core functionality. It also requires `libnuma` (≥v2.0.9) for NUMA-aware memory management and `libtcmalloc_minimal` from Google Performance Tools (gperftools ≥v2.1) for high-performance memory allocation.

## Replicating Our Benchmarks
Get the reordering code:
```
> git clone https://github.com/PASSIONLab/reordering-spgemm.git
> export REORDERING_ROOT=$(pwd)/reordering-spgemm

> cd $REORDERING_ROOT
```
Download the datasets:
```
> export DATA_PATH=$(pwd)/data

> sh $DATA_PATH/download_datasets.sh
> sh $DATA_PATH/unzip_datasets.sh

> mkdir $DATA_PATH/reordering
```

### RCM & ND
Build and install libmtx:
```
> export LIBMTX_HOME=$(pwd)/libmtx
> cd $LIBMTX_HOME

> autoreconf -i -f
> ./configure
> make
> sudo make install
```
Run RCM and ND reordering:
```
> mkdir $DATA_PATH/reordering/rcm_order
> mkdir $DATA_PATH/reordering/nd_order

> export RCM_DATA_PATH=$DATA_PATH/reordering/rcm_order
> export ND_DATA_PATH=$DATA_PATH/reordering/nd_order

> sh reorder_rcm.sh
> sh reorder_nd.sh
```
### AMD, Gray, Degree & SlashBurn
```
> export SPARSEBASE_HOME=$(pwd)/SparseBase
> cd $SPARSEBASE_HOME
```
Build SparseBase:
```
> mkdir build && cd build
> cmake -D_HEADER_ONLY=ON -DCMAKE_BUILD_TYPE=Release -DUSE_CUDA=OFF ..
> make
```
To build SparseBase with AMD, install SuiteSparse first:
```
> git clone https://github.com/DrTimothyAldenDavis/SuiteSparse.git
> cd SuiteSparse 
> mkdir -p build && cd build
> cmake ..
> cmake --build .
> sudo cmake --install .
```
Build SparseBase with AMD:
```
> cd $SPARSEBASE_HOME
> mkdir build && cd build
> cmake .. -D_HEADER_ONLY=ON -DCMAKE_BUILD_TYPE=Release -DUSE_CUDA=OFF -DUSE_AMD_ORDER=ON -DAMD_LIB_DIR=/usr/local/lib -DAMD_INC_DIR=/usr/local/include/suitesparse
> make
```
Run AMD, Gray, Degree, and SlashBurn reordering:
```
> mkdir $DATA_PATH/reordering/amd_order
> export AMD_DATA_PATH=$DATA_PATH/reordering/amd_order

> mkdir $DATA_PATH/reordering/gray_order
> export GRAY_DATA_PATH=$DATA_PATH/reordering/gray_order

> mkdir $DATA_PATH/reordering/degree_order
> export DEGREE_DATA_PATH=$DATA_PATH/reordering/degree_order

> mkdir $DATA_PATH/reordering/slashburn_order
> export SLASHBURN_DATA_PATH=$DATA_PATH/reordering/slashburn_order

> sh "$SPARSEBASE_HOME/scripts/reorder_amd.sh"
> sh "$SPARSEBASE_HOME/scripts/reorder_gray.sh"
> sh "$SPARSEBASE_HOME/scripts/reorder_degree.sh"
> sh "$SPARSEBASE_HOME/scripts/reorder_slashburn.sh"
```
### GP
Build and install METIS:
```
> export METIS_HOME=$(pwd)/Matrix-Partitioning-Utility/METIS
> cd $METIS_HOME

> wget https://karypis.github.io/glaros/files/sw/metis/metis-5.1.0.tar.gz
> tar -xvzf metis-5.1.0.tar.gz

> g++ -g mmio.c main.cpp -o run -std=c++11 -I"${METIS_HOME}/metis-5.1.0/include" -L"${METIS_HOME}/metis-5.1.0/lib" -Wl,-rpath="${METIS_HOME}/metis-5.1.0/lib" -lmetis
```
Run GP reordering:
```
> mkdir $DATA_PATH/reordering/gp_order
> export GP_DATA_PATH=$DATA_PATH/reordering/gp_order

> sh gp_partition.sh
```
Convert GP partition to ordering:
```
> export PARTITION_HOME=$(pwd)/Matrix-Partitioning-Utility
> cd $PARTITION_HOME/Converter

> make
> sh gp_convert.sh
```
### HP
Build and install PaToH:
```
> export PATOH_HOME=$(pwd)/Matrix-Partitioning-Utility/PaToH
> cd $PATOH_HOME

> gcc main.c libpatoh_linux.a -lm
```
Run HP reordering:
```
> mkdir $DATA_PATH/reordering/hp_order
> export HP_DATA_PATH=$DATA_PATH/reordering/hp_order

> sh hp_partition.sh
```
Convert HP partition to ordering:
```
> export PARTITION_HOME=$(pwd)/Matrix-Partitioning-Utility
> cd $PARTITION_HOME/Converter

> make
> sh hp_convert.sh
```
### Rabbit
Build rabbit_order:
```
> export RABBIT_HOME=$(pwd)/rabbit_order
> cd $RABBIT_HOME/demo
> make
```
Run Rabbit reordering:
```
> mkdir $DATA_PATH/reordering/rabbit_order
> export RABBIT_DATA_PATH=$DATA_PATH/reordering/rabbit_order

> sh reorder_rabbit.sh
```
## Citation and Acknowledgments
If you find this code useful, please cite our paper:
```
@article{cluster-wise-spgemm,
  title={Improving SpGEMM Performance Through Matrix Reordering and Cluster-wise Computation},
  author={Islam, Abdullah Al Raqibul and Xu, Helen and Dai, Dong and Bulu{\c{c}}, Ayd{\i}n},
  journal={arXiv preprint arXiv:2507.21253},
  year={2025}
}
```

## References
1. W.-H. Liu and A. H. Sherman. *Comparative analysis of the Cuthill–McKee and the reverse Cuthill–McKee ordering algorithms for sparse matrices.*  
   SIAM Journal on Numerical Analysis, 13(2):198–213, 1976.

2. P. R. Amestoy, T. A. Davis, and I. S. Duff. *Algorithm 837: AMD, an approximate minimum degree ordering algorithm.*  
   ACM Transactions on Mathematical Software (TOMS), 30(3):381–388, 2004.

3. A. George. *Nested dissection of a regular finite element mesh.*  
   SIAM Journal on Numerical Analysis, 10(2):345–363, 1973.

4. G. Karypis and V. Kumar. *A Fast and High Quality Multilevel Scheme for Partitioning Irregular Graphs.*  
   SIAM Journal on Scientific Computing, 20(1):359–392, 1998.  
   doi: [10.1137/S1064827595287997](https://doi.org/10.1137/S1064827595287997)

5. U. V. Catalyurek and C. Aykanat. *Hypergraph-Partitioning-Based Decomposition for Parallel Sparse-Matrix Vector Multiplication.*  
   IEEE Transactions on Parallel and Distributed Systems, 10(7):673–693, 1999.  
   doi: [10.1109/71.780863](https://doi.org/10.1109/71.780863)

6. H. Zhao, T. Xia, C. Li, W. Zhao, N. Zheng, and P. Ren. *Exploring better speculation and data locality in sparse matrix-vector multiplication on Intel Xeon.*  
   2020 IEEE 38th International Conference on Computer Design (ICCD), pp. 601–609, 2020. IEEE.

7. J. Arai, H. Shiokawa, T. Yamamuro, M. Onizuka, and S. Iwamura. *Rabbit Order: Just-in-Time Parallel Reordering for Fast Graph Analysis.*  
   Proceedings of the 2016 IEEE International Parallel and Distributed Processing Symposium (IPDPS), pp. 22–31, 2016. IEEE.  
   doi: [10.1109/IPDPS.2016.15](https://doi.org/10.1109/IPDPS.2016.15)

8. Y. Lim, U. Kang, and C. Faloutsos. *SlashBurn: Graph Compression and Mining beyond Caveman Communities.*  
   IEEE Transactions on Knowledge and Data Engineering, 26(12):3077–3089, 2014.  
   doi: [10.1109/TKDE.2014.2320716](https://doi.org/10.1109/TKDE.2014.2320716)
