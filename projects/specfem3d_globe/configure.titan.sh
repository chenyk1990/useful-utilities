#!/bin/bash

#see: ~/.bashrc
source /etc/profile
module list

mpif90=ftn
mpicc=cc
f90=ftn
cc=cc

## gnu compilers
warn="-Wunused -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow"
flags="-O3 -mcmodel=medium $warn"
cflags=""

##################################################
# with asdf, adios and cuda5 support
# 1. Make sure you download and compiled the asdf library(https://github.com/SeismicData/asdf-library).
# 2. load adios and cuda library from titan system module

./configure --with-asdf ASDF_LIBS="/ccs/home/chenyk/asdf-library/build/lib/libasdf.a" --with-adios --with-cuda=cuda5 --host=x86_64-unknown-linux-gnu MPIF90=$mpif90 F90=$f90 CC=$cc FLAGS_CHECK="$flags" FCFLAGS="" CFLAGS="$cflags" CUDA_INC="$CUDATOOLKIT_HOME/include" CUDA_LIB="$CUDATOOLKIT_HOME/lib64" MPI_INC="$CRAY_MPICH2_DIR/include"

##################################################
# with adios and cuda
# load adios and cuda library from titan system module

#./configure --with-adios --with-cuda=cuda5 --host=x86_64-unknown-linux-gnu MPIF90=$mpif90 F90=$f90 CC=$cc FLAGS_CHECK="$flags" FCFLAGS="" CFLAGS="$cflags" CUDA_INC="$CUDATOOLKIT_HOME/include" CUDA_LIB="$CUDATOOLKIT_HOME/lib64" MPI_INC="$CRAY_MPICH2_DIR/include"

##
## setup
##
echo
echo "modifying mesh_constants_cuda.h..."
sed -i "/ENABLE_VERY_SLOW_ERROR_CHECKING/ c\#undef ENABLE_VERY_SLOW_ERROR_CHECKING" src/gpu/mesh_constants_cuda.h

echo
echo "done"
echo

