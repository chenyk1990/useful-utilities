#!/bin/bash


gfortran getCN1maps.f -o getCN1maps
gfortran getCN1point.f -o getCN1point
gfortran getCN1xyz.f -o getCN1xyz
gfortran getCN1point_modified.f -o getCN1point_modified

rm crust1.xyz
# extract 1D profile at a single point
# lat="$1"
# lon="$2"
# extract a map
#lat1=-90
#lon1=-180

#lat2=90
#lon2=180

lat1=-90
lat2=90
lon1=-180
lon2=180

# interval
ilat=2
ilon=2
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
echo $lon $lat $topo $moho | awk '{printf "%8.2f %8.2f %5.2f %5.2f\n",$1,$2,$3,$4}' >> crust1_whole_d2.txt
done
done

