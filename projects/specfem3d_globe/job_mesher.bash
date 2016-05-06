#!/bin/bash

##PBS -A GEO111
#PBS -A CSC190SPECFEM
#PBS -N SPECFEM3D_mesher
#PBS -j oe
#PBS -o job_mesher.$PBS_JOBID.o

###################################################

## USER PARAMETERS
## chester: gpu compute nodes have 1 GPU card (K20x) and 16-core (interlagos) CPU

#PBS -l walltime=1:00:00
#PBS -l nodes=6

###################################################

cd $PBS_O_WORKDIR

# read DATA/Par_file to get information about the run
# compute total number of processors needed
NPROC_XI=`grep NPROC_XI DATA/Par_file | cut -d = -f 2 `
NPROC_ETA=`grep NPROC_ETA DATA/Par_file | cut -d = -f 2`
NCHUNKS=`grep NCHUNKS DATA/Par_file | cut -d = -f 2 `

# total number of nodes is the product of the values read
numproc=$(( $NCHUNKS * $NPROC_XI * $NPROC_ETA ))

echo "running simulation: `date`"
echo "working directory: `pwd`"

# cleanup
rm -rf OUTPUT_FILES/*
rm -rf DATABASES_MPI/*

# obtain job information
cat $PBS_NODEFILE > OUTPUT_FILES/compute_nodes
echo "$PBS_JOBID" > OUTPUT_FILES/jobid

# stores setup
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

# runs mesher
echo
echo "running mesher on $numproc processors..."
echo "Job start: `date`"

aprun -n $numproc ./bin/xmeshfem3D

echo "Job done: `date`"
echo
echo "see results in directory: OUTPUT_FILES/"
echo

