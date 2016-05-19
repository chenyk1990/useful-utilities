#!/bin/bash

rundir=./

Nrun=3
nodes=$((96*$Nrun)) 
sed -i "s/#PBS -l nodes=.*/#PBS -l nodes=$nodes/" job_solver_orig.bash
sed -i "s/Nrun=.*/Nrun=$Nrun/" job_solver_orig.bash
sed -i "s/nodes=.*/nodes=$nodes/" job_solver_orig.bash

stationlist=`ls /lustre/atlas/proj-shared/geo111/rawdata/stations/ | head -$Nrun`
eventlist=`ls /lustre/atlas/proj-shared/geo111/rawdata/cmt/cmt.GCMT/ |head -$Nrun`

sed -i "s/.*NUMBER_OF_SIMULTANEOUS_RUNS     =.*/NUMBER_OF_SIMULTANEOUS_RUNS     =$Nrun/"  DATA/Par_file
sed -i 's:.*BROADCAST_SAME_MESH_AND_MODEL   =.*:BROADCAST_SAME_MESH_AND_MODEL   = .true.:' DATA/Par_file

id=0
for ievent in ${eventlist[@]}
do
    id=$(($id+1))
    rundir=run`printf "%04i\n" $id`
    mkdir -p $rundir
    mkdir -p $rundir/DATA
    mkdir -p $rundir/DATABASES_MPI
    mkdir -p $rundir/OUTPUT_FILES
    cp -f /lustre/atlas/proj-shared/geo111/rawdata/cmt/cmt.GCMT/$ievent $rundir/DATA/CMTSOLUTION 
done

id=0
for istation in ${stationlist[@]}
do
    id=$(($id+1))
    rundir=run`printf "%04i\n" $id`
    cp -f /lustre/atlas/proj-shared/geo111/rawdata/stations/$istation $rundir/DATA/STATIONS
done

