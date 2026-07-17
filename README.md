# SpGEMM-Baselines
Other six SpGEMM libraries that compared with LeSpGEMM.

## 0. Preparation

```bash
export SPGEMM_BASELINES=/path/to/SpGEMM-Baselines
export SPARSEOPS_ROOT=/path/to/SparseOps
export DATASETS=$SPARSEOPS_ROOT/data
```
## 1.  HaSpGEMM
1. Installation and quick test:
```bash
cd $SPGEMM_BASELINES/HASpGEMM
make

./haspgemm12 $DATASETS/cant/cant.mtx $DATASETS/cant/cant.mtx
```
2. Run scripts：
```bash
DATA_PATH=$DATASETS  sh run_benchmark.sh
```
The results will be written to `haspgemm_results.txt`.

## 2.  AOCL-Sparse
1. Requirement libraries :
- aocl-libs : `amd-blis`、`amd-libflame`、`amd-utils`
These three libraries should be installed in the same directory, such as `/usr/local/lib`.

2. Installation and quick test:
```bash
cd $SPGEMM_BASELINES/aocl-sparse-5.2.2
mkdir build 

cmake -S $SPGEMM_BASELINES/aocl-sparse-5.2.2 -B $SPGEMM_BASELINES/aocl-sparse-5.2.2/build \
  -DCMAKE_AOCL_ROOT=/path/to/requirement/aocl-libs/

cmake --build $SPGEMM_BASELINES/aocl-sparse-5.2.2/build -j"$(nproc)"
```
3. Run scripts：
```bash
MTX_ROOT=$DATASETS sh run_aoclspgemm.sh
```
The results will be written to `aoclspgemm_results.txt`.

## 3.  Kokkos-Kernels
0. Pre-install Kokkos for requirement:
```bash
cd $SPGEMM_BASELINES/kokkos

cmake -B builddir -S . \
-DCMAKE_INSTALL_PREFIX=$SPGEMM_BASELINES/kokkos/install \
-DKokkos_ENABLE_SERIAL=On \
-DKokkos_ENABLE_OPENMP=On

cmake --build builddir -j

cmake --install builddir
```
1. Compile and install kokkos-kernels:
```bash
cd $SPGEMM_BASELINES/kokkos-kernels

cmake -S . -B build \
-DCMAKE_BUILD_TYPE=Release \
-DKokkos_DIR=$SPGEMM_BASELINES/kokkos/install/lib64/cmake/Kokkos \
-DKokkosKernels_INST_DOUBLE=On \
-DKokkosKernels_ENABLE_PERFTESTS=ON \
-DKokkosKernels_ENABLE_TESTS_AND_PERFSUITE=ON \
-DKokkosKernels_ENABLE_TESTS=ON

cmake --build build -j
make -j
```
Quick test for SpGEMM:
```bash
./perf_test/sparse/sparse_spgemm --amtx $DATASETS/cant/cant.mtx --repeat 5 --openmp 128
```
2. Run scripts for 20 matrices:
```bash
./run_spgemm_benchmark.sh -d casedataset.txt -p $DATASETS -r 5 -t 128 -o Kokkos_results.txt
```
The results will be written to `Kokkos_results.txt`.

## 4.  CSeg
1. Installation and run scripts for testing:
```bash
cd $SPGEMM_BASELINES/CSeg
make
DATA_PATH=$DATASETS sh run_benchmark.sh
```
The results will be written to `CSeg_results.txt`.

## 5. SuiteSparse: GraphBLAS
1. Compile and Installation:
```bash
cd $SPGEMM_BASELINES/GraphBLAS
mkdir build && cd build
cmake ..
make -j
cd .. && make demos 
```
2. Quick test:
```
./build/spgemm_demo $DATASETS/cant/cant.mtx $DATASETS/cant/cant.mtx 
```
3. Run scripts:
```bash
DATA_PATH=$DATASETS sh run_benchmark.sh
```
The results will be written to `GraphBLAS_spgemm_results.txt`.

## 6. Hierarchical $AA^T$ SpGEMM
1. Compilation and Installation:
```bash
cd $SPGEMM_BASELINES/SpGEMM-exp/clusterwise-spgemm-main

make clean && make sample_hw

export CLUSTERWISE_ROOT=$SPGEMM_BASELINES/SpGEMM-exp/clusterwise-spgemm-main
export CLOSE_PAIR_DATA_PATH=$SPARSEOPS_ROOT/data/reordering/close_pair
```
2. Quick test:
```bash
./bin/HierarchicalClusterSpGEMM_hw text $DATASETS/cant/cant.mtx $DATASETS/cant/cant.mtx $CLOSE_PAIR_DATA_PATH/cant.mtx
```
2. Run scripts:
```bash
DATA_PATH=$DATASETS sh scripts/h_clusterwise/2_hierarchical_spgemm_from_file.sh
```
The results will be printed on screen or you can redirect the output by:
```bash
DATA_PATH=$DATASETS sh scripts/h_clusterwise/2_hierarchical_spgemm_from_file.sh > Hierarchical_SpGEMM_Results.txt
```

