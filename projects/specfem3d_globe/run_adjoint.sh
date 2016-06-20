#!/bin/bash

rundir=./

Nrun=10
nodes=$((96*$Nrun)) 
sed -i "s/#PBS -l nodes=.*/#PBS -l nodes=$nodes/" job_solver_multishots.bash
sed -i "s/Nrun=.*/Nrun=$Nrun/" job_solver_multishots.bash
sed -i "s/nodes=.*/nodes=$nodes/" job_solver_multishots.bash

#eventlist=`ls /lustre/atlas/proj-shared/geo111/rawdata/asdf/raw_obsd/| head -$(($Nrun)) |cut -c-8 | tail -$Nrun`
eventlist=`ls /lustre/atlas/proj-shared/csc190/specfem/chenyk/asdf_raw | head -$(($Nrun+1)) |cut -c-8|sed  "5d"`

sed -i "s/.*SIMULATION_TYPE                 =.*#/SIMULATION_TYPE                 = 3 #/" DATA/Par_file
sed -i "s/.*SAVE_FORWARD                    =.*#/SAVE_FORWARD                    =.False. #/" DATA/Par_file
sed -i "s/.*NUMBER_OF_SIMULTANEOUS_RUNS     =.*/NUMBER_OF_SIMULTANEOUS_RUNS     =$Nrun/"  DATA/Par_file
sed -i 's:.*BROADCAST_SAME_MESH_AND_MODEL   =.*:BROADCAST_SAME_MESH_AND_MODEL   = .true.:' DATA/Par_file
sed -i "s/.*READ_ADJSRC_ASDF                =.*/READ_ADJSRC_ASDF                =.true./" DATA/Par_file

id=0
for ievent in ${eventlist[@]}
do
    id=$(($id+1))
    rundir=run`printf "%04i\n" $id`
    mkdir -p $rundir/SEM 
    cp -f /lustre/atlas/proj-shared/csc190/specfem/chenyk/princeton/processing/proc/adjsum/$ievent.adjoint.h5 $rundir/SEM/adjoint.h5 
    pypaw-generate_stations_asdf /lustre/atlas/proj-shared/csc190/specfem/chenyk/princeton/processing/proc/adjsum/$ievent.adjoint.h5 
    mv STATIONS_ADJOINT $rundir/DATA/
done
