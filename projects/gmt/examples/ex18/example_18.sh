#!/bin/bash
#		GMT EXAMPLE 18
#		$Id: example_18.sh 15178 2015-11-06 10:45:03Z fwobbe $
#
# Purpose:	Illustrates volumes of grids inside contours and spatial
#		selection of data
# GMT progs:	gmtset, gmtselect, gmtspatial, grdclip, grdcontour, grdgradient, grdimage
# GMT progs:	grdmath, grdvolume, makecpt, pscoast, psscale, pstext, psxy
# Unix progs:	$AWK, cat, rm
#
ps=example_18.ps

# Use spherical gmt projection since SS data define on sphere
gmt gmtset PROJ_ELLIPSOID Sphere FORMAT_FLOAT_OUT %g

# Define location of Pratt seamount and the 400 km diameter
echo "-142.65 56.25 400" > pratt.d

# First generate gravity image w/ shading, label Pratt, and draw a circle
# of radius = 200 km centered on Pratt.

gmt makecpt -Crainbow -T-60/60/120 -Z > grav.cpt
gmt grdgradient AK_gulf_grav.nc -Nt1 -A45 -GAK_gulf_grav_i.nc
gmt grdimage AK_gulf_grav.nc -IAK_gulf_grav_i.nc -JM5.5i -Cgrav.cpt -B2f1 -P -K -X1.5i \
	-Y5.85i > $ps
gmt pscoast -RAK_gulf_grav.nc -J -O -K -Di -Ggray -Wthinnest >> $ps
gmt psscale -DJBC+o0/0.4i+w4i/0.15i+h -R -J -Cgrav.cpt -Bx20f10 -By+l"mGal" -O -K >> $ps
$AWK '{print $1, $2, "Pratt"}' pratt.d | gmt pstext -R -J -O -K -D0.1i/0.1i \
	-F+f12p,Helvetica-Bold+jLB >> $ps
gmt psxy pratt.d -R -J -O -K -SE- -Wthinnest >> $ps

# Then draw 10 mGal contours and overlay 50 mGal contour in green

gmt grdcontour AK_gulf_grav.nc -J -C20 -B2f1 -BWSEn -O -K -Y-4.85i >> $ps
# Save 50 mGal contours to individual files, then plot them
gmt grdcontour AK_gulf_grav.nc -C10 -L49/51 -Dsm_%d_%c.txt
gmt psxy -R -J -O -K -Wthin,green sm_*.txt >> $ps
gmt pscoast -R -J -O -K -Di -Ggray -Wthinnest >> $ps
gmt psxy pratt.d -R -J -O -K -SE- -Wthinnest >> $ps
rm -f sm_*_O.txt	# Only consider the closed contours

# Now determine centers of each enclosed seamount > 50 mGal but only plot
# the ones within 200 km of Pratt seamount.

# First determine mean location of each closed contour and
# add it to the file centers.d

gmt gmtspatial -Q -fg sm_*_C.txt > centers.d

# Only plot the ones within 200 km

gmt gmtselect -C200k/pratt.d centers.d -fg | gmt psxy -R -J -O -K -SC0.04i -Gred -Wthinnest >> $ps
gmt psxy -R -J -O -K -ST0.1i -Gyellow -Wthinnest pratt.d >> $ps

# Then report the volume and area of these seamounts only
# by masking out data outside the 200 km-radius circle
# and then evaluate area/volume for the 50 mGal contour

gmt grdmath -R `$AWK '{print $1, $2}' pratt.d` SDIST = mask.nc
gmt grdclip mask.nc -Sa200/NaN -Sb200/1 -Gmask.nc
gmt grdmath AK_gulf_grav.nc mask.nc MUL = tmp.nc
area=`gmt grdvolume tmp.nc -C50 -Sk | cut -f2`
volume=`gmt grdvolume tmp.nc -C50 -Sk | cut -f3`

gmt psxy -R -J -A -O -K -L -Wthin -Gwhite >> $ps << END
-148.5	52.75
-141	52.75
-141	53.75
-148.5	53.75
END
gmt pstext -R -J -O -F+f14p,Helvetica-Bold+jLM >> $ps << END
-148 53.08 Areas: $area km@+2@+
-148 53.42 Volumes: $volume mGal\264km@+2@+
END

# Clean up

rm -f grav.cpt sm_*.txt *_i.nc tmp.nc mask.nc pratt.d center* gmt.conf
