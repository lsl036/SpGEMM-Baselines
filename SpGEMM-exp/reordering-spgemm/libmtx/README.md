# Libmtx
Libmtx is a C library and collection of utility programs for working with objects in the Matrix Market file
format, including dense and sparse matrices and vectors.

  Copyright (C) 2022 James D. Trotter

  Copying and distribution of this file, with or without modification,
  are permitted in any medium without royalty provided the copyright
  notice and this notice are preserved.

See the file INSTALL for instructions on how to build and install. Libmtx is free software. See the file COPYING for copying conditions.

## Build and Install
We used Libmtx to perform matrix reordering with `Reverse Cuthill-McKee (RCM)` and `Nested Dissection (ND)`. Use this command to build and install libmtx.

```
> autoreconf -i -f
> ./configure
> make
> sudo make install
```

After this, libmtx will be installed in `/usr/local/lib/`.

## Sample Command
We used `mtxreorder` command to compute RCM and ND orderings. The following commands generate the ordering files from a matrix file in Matrix Market format.
```
> mtxreorder --verbose /path-to-file/dataset.mtx --ordering=rcm --rowperm-path=/path-to-file/dataset.rcmorder -q
> mtxreorder --verbose /path-to-file/dataset.mtx --ordering=nd --rowperm-path=/path-to-file/dataset.ndorder -q
```