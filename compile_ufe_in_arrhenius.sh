#!/bin/bash

# Load required modules (adapt as needed)
# I have not checked the modules! job to be done, compare with Tim's code
module purge
module load buildenv-gcc/2025b-eb
module load netCDF/4.9.3 #netCDF-HDF5/4.9.2-1.12.2-hpc1
module load METIS/5.1.0 
module load CMake/3.31.8 Ninja/1.13.0 HDF5/1.14.6
# this modules I are needed and being loaded by buildenv-gcc/2025b-eb
#OpenBLAS/0.3.30 FlexiBLAS/3.4.5 FFTW/3.3.10 FFTW.MPI/3.3.10 ScaLAPACK/2.2.2-fb

# PETSc paths
export PETSC_DIR=/nobackup/proj/disk/bolinc/personal/frare/petsc/petsc_install/3.22.4_gcc2025b-eb
export PETSC_ARCH=  # Leave unset if you're using a flat install (recommended)
export PATH=$PETSC_DIR/bin:$PATH
export LD_LIBRARY_PATH=$PETSC_DIR/lib:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=$PETSC_DIR/lib/pkgconfig:$PKG_CONFIG_PATH

# Set NetCDF directory from the module information
export NETCDF_DIR=/software/sse2/el9_epyc9005/easybuild/prefix/software/netCDF/4.9.3-gompi-2025b
export NETCDF_INC=${NETCDF_DIR}/include
export NETCDF_LIB=${NETCDF_DIR}/lib

BUILD_TYPE=${1:-dev}        # default: dev
BUILD_MODE=${2:-clean}      # default: clean

echo ""
echo ">>> Compiling UFEMISM ($BUILD_TYPE, $BUILD_MODE)"
echo ""

cd UFEMISM2.0_porting/

# Clean or incremental build handling
if [ "$BUILD_MODE" == "clean" ]; then
    echo ">> Performing clean build"
    rm -rf build
elif [ "$BUILD_MODE" == "changed" ]; then
    echo ">> Performing incremental build"
    rm -f build/CMakeCache.txt
else
    echo "Unknown build mode: $BUILD_MODE"
    exit 1
fi

# Create build directory if missing
mkdir -p build

COMMON_FLAGS="-Wall;-ffree-line-length-none;-cpp;-fimplicit-none;-g;-I${NETCDF_INC}"

if [ "$BUILD_TYPE" == "dev" ]; then
    cmake -B build -S . -G Ninja \
      -DBUILD_UFEMISM=ON \
      -DPETSC_DIR=$PETSC_DIR \
      -DDO_ASSERTIONS=ON \
      -DDO_RESOURCE_TRACKING=ON \
      -DEXTRA_Fortran_FLAGS="${COMMON_FLAGS};-fdiagnostics-color=always;-Og;-Werror=implicit-interface;-fcheck=all;-fbacktrace;-finit-real=nan;-finit-integer=-42;-finit-character=33" \
      -DCMAKE_EXE_LINKER_FLAGS="-L${PETSC_DIR}/lib -lpetsc -L${NETCDF_LIB} -lnetcdff -lnetcdf"
elif [ "$BUILD_TYPE" == "perf" ]; then
    cmake -B build -S . -G Ninja \
      -DBUILD_UFEMISM=ON \
      -DPETSC_DIR=$PETSC_DIR \
      -DDO_ASSERTIONS=OFF \
      -DDO_RESOURCE_TRACKING=OFF \
      -DEXTRA_Fortran_FLAGS="$COMMON_FLAGS -O3"
else
    echo "Unknown build type: $BUILD_TYPE"
    exit 1
fi

# Build the project using Ninja
ninja -C build -v

# Copy binary to top directory
BINNAME=UFEMISM_program
BINPATH=build/src/UFEMISM/${BINNAME}

if [ "$BUILD_TYPE" == "dev" ]; then
    cp ${BINPATH} ${BINNAME}_dev
    cp ${BINNAME}_dev ${BINNAME}
elif [ "$BUILD_TYPE" == "perf" ]; then
    cp ${BINPATH} ${BINNAME}_perf
    cp ${BINNAME}_perf ${BINNAME}
fi

echo ""
echo "Compilation finished. Binary: $BINNAME"