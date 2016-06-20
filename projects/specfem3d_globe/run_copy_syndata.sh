#!/bin/bash

Nrun=10


#eventlist=`ls /lustre/atlas/proj-shared/geo111/rawdata/asdf/raw_obsd/| head -$(($Nrun)) | cut -c-8`
eventlist=`ls /lustre/atlas/proj-shared/csc190/specfem/chenyk/asdf_raw | head -$(($Nrun+1)) |cut -c-8 | sed "5d"`

id=0
for ievent in ${eventlist[@]}
do
      id=$(($id+1))
      rundir=run`printf "%04i\n" $id`
      cp $rundir/OUTPUT_FILES/synthetic.h5 /lustre/atlas/proj-shared/csc190/specfem/chenyk/princeton/processing/asdf_syn/$ievent.synthetic.h5
done
