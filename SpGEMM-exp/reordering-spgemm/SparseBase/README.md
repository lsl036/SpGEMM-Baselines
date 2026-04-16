# SparseBase

SparseBase is a library built in C++ that encapsulates, preprocesses, and performs I/O operations on sparse data structures seamlessly and optimally to provide a backbone for algorithms that use these structures.

It is designed with HPC (High Performance Computing) usage at the forefront. It is meant as a container of sparse objects such as tensors, graphs, multi-graphs, and hypergraphs. It mainly focuses on re-ordering, partitioning, and coarsening sparse objects. Also, it supports many different formats (representations) of sparse data and enables smooth conversion between formats.

The library is still in early stages of development. As a result, the API is not stable and likely to change in the near future.

:rocket: [Installation & Basic Usage](https://sparcityeu.github.io/sparsebase/pages/getting_started.html)

:computer: [Source Code](https://github.com/sparcityeu/sparsebase)

:books: [Documentation](https://sparcityeu.github.io/SparseBase/)

:scroll: [License](https://sparcityeu.github.io/sparsebase/pages/license.html)

:heart: [Contribution Guide](https://sparcityeu.github.io/sparsebase/pages/contributing/index.html)

## Build SparseBase Code
```
> mkdir build && cd build
> cmake -D_HEADER_ONLY=ON -DCMAKE_BUILD_TYPE=Release -DUSE_CUDA=OFF ..
> make
```
## Build SparseBase With AMD
### Install SuiteSparse
```
> mkdir -p build && cd build
> cmake ..
> cmake --build .
> sudo cmake --install .
or
> sudo make install
```
### Install SuiteSparse only in SuiteSparse/lib (not /usr/local/lib)
```
> mkdir -p build && cd build
> cmake -DCMAKE_INSTALL_PREFIX=.. ..
> cmake --build .
> cmake --install .
```
## Build SparseBase
### Build SparseBase Code
```
> mkdir build && cd build
> cmake .. -D_HEADER_ONLY=ON -DCMAKE_BUILD_TYPE=Release -DUSE_CUDA=OFF -DUSE_AMD_ORDER=ON -DAMD_LIB_DIR=/usr/local/lib -DAMD_INC_DIR=/usr/local/include/suitesparse
> make
```

### Build SparseBase Code: SuiteSparse only in SuiteSparse/lib (not /usr/local/lib)
```
> mkdir build && cd build
> cmake .. -D_HEADER_ONLY=ON -DCMAKE_BUILD_TYPE=Release -DUSE_CUDA=OFF -DUSE_AMD_ORDER=ON -DAMD_LIB_DIR=/global/homes/r/raqib/SuiteSparse/lib64 -DAMD_INC_DIR=/global/homes/r/raqib/SuiteSparse/include/suitesparse
> make
```

## Run Examples
```
> ./examples/gray_order/gray_order "{DATA_PATH}/2cubes_sphere/2cubes_sphere.mtx" "{DATA_PATH}/2cubes_sphere/2cubes_sphere.grayorder"
> ./examples/amd_order/amd_order "{DATA_PATH}/2cubes_sphere/2cubes_sphere.mtx" "{DATA_PATH}/2cubes_sphere/2cubes_sphere.amdorder"
```