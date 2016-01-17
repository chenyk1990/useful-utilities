
## submit jobs
sbatch dblend.job
##
checking status
squeue -u chenyk
showq  -u chenyk
## apply node
idev -A TCCS -N 2 -n 32 -t 1:00:00 -p normal
ibrun -run MPI program

## denotions
N: number of nodes
n: number of MPI task
A: project name

When using pscons
key point:
set RSF_THREADS=16*N
set RSF_CLUSTER=`host.sh`
pscons then works

