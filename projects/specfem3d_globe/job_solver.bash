#!/bin/bash

##PBS -A GEO111
#PBS -A CSC190SPECFEM
#PBS -N SPECFEM3D_solver
#PBS -j oe
#PBS -o job_solver.$PBS_JOBID.o

###################################################

## USER PARAMETERS
## chester: gpu compute nodes have 1 GPU card (K20x) and 16-core (interlagos) CPU

#PBS -l walltime=1:00:00
#PBS -l nodes=96

###################################################

cd $PBS_O_WORKDIR

NPROC_XI=`grep NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep NPROC_ETA DATA/Par_file | cut -d = -f 2`
NCHUNKS=`grep NCHUNKS DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numproc=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))

echo "running simulation: `date`"
echo "directory: `pwd`"
echo

# obtain job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# runs simulation
echo
echo "running solver..."
echo `date`
aprun -n $numproc -N1 ./bin/xspecfem3D

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done: `date`"

