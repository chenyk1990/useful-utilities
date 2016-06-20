#!/bin/bash

Nrun=10

#ls /lustre/atlas/proj-shared/geo111/rawdata/asdf/raw_obsd/| head -$(($Nrun)) |cut -c-8 >cmtlist
ls /lustre/atlas/proj-shared/geo111/Wenjie/DATA_EBRU/raw/| head -$(($Nrun)) |cut -c-8 >cmtlist
#eventlist=`ls /lustre/atlas/proj-shared/geo111/Wenjie/DATA_EBRU/raw/| head -$(($Nrun)) |cut -c-8 `
eventlist=`ls /lustre/atlas/proj-shared/csc190/specfem/chenyk/asdf_raw | head -$(($Nrun+1)) |cut -c-8 | sed "5d"`
ls /lustre/atlas/proj-shared/csc190/specfem/chenyk/asdf_raw | head -$(($Nrun+1)) |cut -c-8 | sed "5d" >cmtlist
#eventlist=`ls /lustre/atlas/proj-shared/geo111/rawdata/asdf/raw_obsd/| head -$(($Nrun))|cut -c-8 `

for ievent in ${eventlist[@]}
do      
#cp /lustre/atlas/proj-shared/geo111/rawdata/asdf/raw_obsd/$ievent.raw_observed.h5 /lustre/atlas/proj-shared/csc190/specfem/chenyk/princeton/processing/asdf_raw/$ievent.observed.h5  
cp  /lustre/atlas/proj-shared/geo111/Wenjie/DATA_EBRU/raw/$ievent.observed.h5 /lustre/atlas/proj-shared/csc190/specfem/chenyk/princeton/processing/asdf_raw/
done
