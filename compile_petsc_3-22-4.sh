#!/bin/bash

# This is a compiled version of PETSc 3.22.4
# for use with UFEMISM on Tetralith.
# It is configured for performance with the GCC 2022a compiler.
# It was needed to checkout to get matGetRow_ functionality. (git checkout v3.22.4)
# be sure which version is using: grep "^#define PETSC_VERSION_" include/petscversion.h
# before installing PETSc.

./configure \
  --prefix=/proj/bolinc/users/x_frare/petsc/petsc_install/3.22.4_gcc2022a_opt_g \
  --with-cc=mpicc --with-cxx=mpicxx --with-fc=mpif90 \
  COPTFLAGS="-O2 -g -march=native" \
  CXXOPTFLAGS="-O2 -g -march=native" \
  FOPTFLAGS="-O2 -g -march=native" \
  --with-debugging=0 \
  --with-log=1 \
  --with-blas-lib="-lopenblas" \
  --with-lapack-lib="-lopenblas" \
  --with-scalapack=1 \
  --with-metis=1 \
  --with-fftw=1 \
  --with-hdf5=1 \
  --with-hwloc=1 \
  --with-suitesparse=1

make all
make install

nm -g /proj/bolinc/users/x_frare/petsc/petsc_install/3.22.4_gcc2022a_opt_g/lib/libpetsc.so | grep -i matgetrow_