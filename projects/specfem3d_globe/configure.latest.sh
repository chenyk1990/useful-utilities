#!/bin/bash

#see: ~/.bashrc
#source /etc/profile
# cuda
#module load cudatoolkit # cuda cudatools cudasdk
## cray compilers
#module unload PrgEnv-gnu
#module load PrgEnv-cray
## gnu compilers
#module unload PrgEnv-cray
#module load PrgEnv-gnu
#module load adios/1.5.0
module list

mpif90=ftn
mpicc=cc
f90=ftn
cc=cc

## cray compilers
#flags="-eF -em -rm"
#cflags="-h list=m"
#mpi_inc="-I/opt/cray/mpt/default/xt/gemini/mpich2-cray/73/include"
#mpi_inc="-I/opt/cray/mpt/default/gni/mpich2-cray/47/include/"

## gnu compilers
warn="-Wunused -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow"
flags="-O3 -mcmodel=medium $warn"
cflags=""
#flags="-std=gnu -fimplicit-none -frange-check -O0 -g -fbacktrace -fbounds-check -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow " # -mcmodel=medium
        # useful for debugging: add -ffpe-trap=overflow,zero -fbacktrace -fbounds-check
#flags="-O2 -mcmodel=medium -g -fbacktrace"
#mpi_inc="-I/opt/cray/mpt/default/gni/mpich2-gnu/47/include/"
#mpi_inc="-I/opt/cray/mpt/default/gni/mpich2-gnu/49/include/"

#adios_link=`/ccs/home/mpbl/adios-1.5.0-chester/bin/adios_config -lf`
#adios_inc=`/ccs/home/mpbl/adios-1.5.0-chester/bin/adios_config -cf`
#adios_c_inc=`/ccs/home/mpbl/adios-1.5.0-chester/bin/adios_config -c`

adios_link=`adios_config -lf`
adios_inc=`adios_config -cf`
adios_c_inc=`adios_config -c`

#cuda_inc="/opt/nvidia/cudatoolkit/5.5.20-1.0402.7700.8.1/include"
#cuda_lib="/opt/nvidia/cudatoolkit/5.5.20-1.0402.7700.8.1/lib64"

# w/ adios 
#./configure --with-cuda=cuda5 --with-adios --host=x86_64-unknown-linux-gnu MPIF90=$mpif90 F90=$f90 CC=$cc FLAGS_CHECK="$flags" FCFLAGS="" CFLAGS="$cflags" CUDA_INC="$CUDATOOLKIT_HOME/include" CUDA_LIB="$CUDATOOLKIT_HOME/lib64" MPI_INC="$CRAY_MPICH2_DIR/include"

# w/ adios 
#./configure --with-cuda=cuda5 --with-adios --host=x86_64-unknown-linux-gnu MPIF90=$mpif90 F90=$f90 CC=$cc FLAGS_CHECK="$flags" FCFLAGS="" CFLAGS="$cflags" CUDA_INC="$CUDATOOLKIT_HOME/include" CUDA_LIB="$CUDATOOLKIT_HOME/lib64" MPI_INC="$CRAY_MPICH2_DIR/include"

# w/ asdf 
#./configure --with-asdf ASDF_LIBS="./libasdf.a" --with-adios --with-cuda=cuda5 --host=x86_64-unknown-linux-gnu MPIF90=$mpif90 F90=$f90 CC=$cc FLAGS_CHECK="$flags" FCFLAGS="" CFLAGS="$cflags" CUDA_INC="$CUDATOOLKIT_HOME/include" CUDA_LIB="$CUDATOOLKIT_HOME/lib64" MPI_INC="$CRAY_MPICH2_DIR/include" 
./configure --with-asdf ASDF_LIBS="/ccs/home/chenyk/asdf-library/build/lib/libasdf.a" --with-adios --with-cuda=cuda5 --host=x86_64-unknown-linux-gnu MPIF90=$mpif90 F90=$f90 CC=$cc FLAGS_CHECK="$flags" FCFLAGS="" CFLAGS="$cflags" CUDA_INC="$CUDATOOLKIT_HOME/include" CUDA_LIB="$CUDATOOLKIT_HOME/lib64" MPI_INC="$CRAY_MPICH2_DIR/include" 
#./configure --with-adios --with-cuda=cuda5 --host=x86_64-unknown-linux-gnu MPIF90=$mpif90 F90=$f90 CC=$cc FLAGS_CHECK="$flags" FCFLAGS="" CFLAGS="$cflags" CUDA_INC="$CUDATOOLKIT_HOME/include" CUDA_LIB="$CUDATOOLKIT_HOME/lib64" MPI_INC="$CRAY_MPICH2_DIR/include" 

# w/ asdf 
#./configure --with-asdf ASDF_LIBS="/lustre/atlas/proj-shared/geo111/test_asdf/specfem3d_globe/libasdf.a" --with-adios --with-cuda=cuda5 --host=x86_64-unknown-linux-gnu MPIF90=$mpif90 F90=$f90 CC=$cc FLAGS_CHECK="$flags" FCFLAGS="" CFLAGS="$cflags" CUDA_INC="$CUDATOOLKIT_HOME/include" CUDA_LIB="$CUDATOOLKIT_HOME/lib64" MPI_INC="$CRAY_MPICH2_DIR/include"

#./configure -host=x86_64-unknown-linux-gnu MPIFC=$mpif90 MPICC=$mpicc FC=$f90 CC=$cc \
#FLAGS_CHECK="$flags" FLAGS_NO_CHECK="$flags" MPI_INCLUDES=$mpi_inc \
# --with-adios \
#ADIOS_INC="$adios_inc" \
#ADIOS_LIB="$adios_link" \
#CPPFLAGS=$mpi_inc
#--with-cuda=cuda5 \
#CUDA_INC=$cuda_inc \
#CUDA_LIB=$cuda_lib \
#CPPFLAGS=$mpi_inc

##
## setup
##
echo
echo "modifiying constants.h..."

# file to change
file=./setup/constants.h

## ETOPO1
#echo "ETOPO1..."
#sed -i "s:NX_BATHY.*:NX_BATHY = 21600,NY_BATHY = 10800:g" $file
#sed -i "s:RESOLUTION_TOPO_FILE.*:RESOLUTION_TOPO_FILE = 1:g" $file
#sed -i "s:PATHNAME_TOPO_FILE.*:PATHNAME_TOPO_FILE = 'DATA/topo_bathy/ETOPO1.xyz':g" $file
#sed -i "s:SMOOTH_CRUST.*:SMOOTH_CRUST = .false.:g" $file

## kernel outputs (isotropic kernels)
echo "kernel outputs..."
#sed -i "s:ANISOTROPIC_KL.*:ANISOTROPIC_KL = .true.:g" $file
#sed -i "s:SAVE_TRANSVERSE_KL.*:SAVE_TRANSVERSE_KL = .true.:g" $file
#sed -i "s:APPROXIMATE_HESS_KL.*:APPROXIMATE_HESS_KL = .true.:g" $file

echo
echo "modifying mesh_constants_cuda.h..."
sed -i "/ENABLE_VERY_SLOW_ERROR_CHECKING/ c\#undef ENABLE_VERY_SLOW_ERROR_CHECKING" src/gpu/mesh_constants_cuda.h

echo
echo "done"
echo

