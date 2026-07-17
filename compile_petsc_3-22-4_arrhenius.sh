#!/bin/bash

# This is a compiled version of PETSc 3.22.4
# for use with UFEMISM on Arrhenius.
# It is configured for performance with the GCC 2022a compiler.
# It was needed to checkout to get matGetRow_ functionality. (git checkout v3.22.4)
# be sure which version is using: grep "^#define PETSC_VERSION_" include/petscversion.h
# before installing PETSc.

./configure \
  --prefix=/nobackup/proj/disk/bolinc/personal/frare/petsc/petsc_install/3.22.4_gcc2025b-eb \
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
  --with-suitesparse=0

make all
make install

nm -g /nobackup/proj/disk/bolinc/personal/frare/petsc/petsc_install/3.22.4_gcc2025b-eb/lib/libpetsc.so | grep -i matgetrow_
# 0000000000c3a520 T matgetrow_