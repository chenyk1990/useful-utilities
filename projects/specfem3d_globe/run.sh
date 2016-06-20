#!/bin/bash

rundir=./

Nrun=10
nodes=$((96*$Nrun)) 
sed -i "s/#PBS -l nodes=.*/#PBS -l nodes=$nodes/" job_solver_multishots.bash
sed -i "s/Nrun=.*/Nrun=$Nrun/" job_solver_multishots.bash
sed -i "s/nodes=.*/nodes=$nodes/" job_solver_multishots.bash

stationlist=`ls /lustre/atlas/proj-shared/geo111/rawdata/stations/ | head -$Nrun`
#eventlist=`ls /lustre/atlas/proj-shared/geo111/rawdata/cmt/cmt.GCMT/ |head -$Nrun`
eventlist=`ls /lustre/atlas/proj-shared/csc190/specfem/chenyk/asdf_raw | head -$(($Nrun+1)) |cut -c-8|sed  "5d"`

sed -i "s/.*SIMULATION_TYPE                 =.*#/SIMULATION_TYPE                 = 1 #/" DATA/Par_file
sed -i "s/.*SAVE_FORWARD                    =.*#/SAVE_FORWARD                    =.True. #/" DATA/Par_file
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
for ievent in ${eventlist[@]}
do
    id=$(($id+1))
    rundir=run`printf "%04i\n" $id`
    cp -f /lustre/atlas/proj-shared/geo111/rawdata/stations/$ievent.STATIONS $rundir/DATA/STATIONS
done
rundir=run0001
#cp -rf DATABASES_MPI/* $rundir/DATABASES_MPI/
#cp -rf OUTPUT_FILES/* $rundir/OUTPUT_FILES/ 

