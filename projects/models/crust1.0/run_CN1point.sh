#!/bin/bash
# shell script to run getCN1point and orgnize the output format for matlab input# Yunfeng Chen, Oct. 7th, 2014, Global Seismology Group, University of Alberta 

gfortran getCN1maps.f -o getCN1maps
gfortran getCN1point.f -o getCN1point
gfortran getCN1xyz.f -o getCN1xyz
gfortran getCN1point_modified.f -o getCN1point_modified

rm crust1.xyz
# extract 1D profile at a single point
# lat="$1"
# lon="$2"
# extract a map
lat1=48
lon1=-121

lat2=58
#lat2=50
lon2=-108
#lon2=-119

# interval
ilat=0.5
ilon=0.5
currdir="/home/yunfeng/20_30/research/crust1.0"
for lat in $(seq $lat1 $ilat $lat2)
do
for lon in $(seq $lon1 $ilon $lon2)
do
echo "Extracting Lat=$lat Lon=$lon" 
./getCN1point << ! > mod.txt
$lon $lat
q
!
topo=($(grep topography mod.txt | awk '{print $2}'))
#moho=($(more +13 mod.txt | head -1 | awk '{print -$4}')) #Yunfeng's
moho=($(tail -n 3 mod.txt | head -1 | awk '{print -$4}'))
echo $lon $lat $topo $moho | awk '{printf "%8.2f %8.2f %5.2f %5.2f\n",$1,$2,$3,$4}' >> crust1.xyz
done
done


