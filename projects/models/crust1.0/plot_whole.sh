#!/bin/bash

PSNAME="H_crust1.ps"

gmtset ELLIPSOID Clarke-1866
gmtset PAPER_MEDIA A4
gmtset BASEMAP_TYPE FANCY PLOT_DEGREE_FORMAT ddd:mm:ssF GRID_CROSS_SIZE_PRIMARY 0.05i

range="-121/-109/48/58"
#range="-180/180/-90/90"
# awk '{print $1, $2, $3}' crust1.0.model | surface -Gtemp.grd -I0.05/0.05 -R$range -T0 -Lld -Lud
awk '{print $1, $2, $4}' crust1.xyz | surface -Gtemp.grd -I0.05/0.05 -R$range -T0 -Lld -Lud
#grdimage temp.grd -JG90/30/5i -K >$PSNAME
# grdmath -V mask.grd temp.grd MUL = H-map.grd
#grd2cpt temp.grd -Crainbow -I  -E50 -D > my_crust1.cpt
#makecpt -T34/54/0.2 -Crainbow -I -D  > my_crust1.cpt
#psbasemap -JG90/30/5i -R$range -B -K -V -P > $PSNAME

#grdimage temp.grd -Cmy_crust1.cpt -JG90/30/5i -R  -K -O -V >> $PSNAME
evince $PSNAME &

# ps2raster $PSNAME -A -Tef -V 
#gv $PSNAME
