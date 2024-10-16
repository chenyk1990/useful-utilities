REM             GMT EXAMPLE 29
REM             $Id: example_29.bat 15178 2015-11-06 10:45:03Z fwobbe $
REM
REM Purpose:	Illustrates spherical surface gridding with Green's function of splines
REM GMT progs:	makecpt, grdcontour, grdgradient, grdimage, grdmath greenspline, psscale, pstext
REM DOS calls:	del, echo
REM

echo GMT EXAMPLE 29
set ps=example_29.ps

REM This example uses 370 radio occultation data for Mars to grid the topography.
REM Data and information from Smith, D. E., and M. T. Zuber (1996), The shape of
REM Mars and the topographic signature of the hemispheric dichotomy, Science, 271, 184–187.

REM Make Mars PROJ_ELLIPSOID given their three best-fitting axes:
set a=3399.472
set b=3394.329
set c=3376.502
gmt grdmath -Rg -I4 -r X COSD %a% DIV DUP MUL X SIND %b% DIV DUP MUL ADD Y COSD DUP MUL MUL Y SIND %c% DIV DUP MUL ADD SQRT INV = PROJ_ELLIPSOID.nc

REM  Do both Parker and Wessel/Becker solutions (tension = 0.9975)
gmt greenspline -RPROJ_ELLIPSOID.nc mars370.in -D4 -Sp -Gmars.nc
gmt greenspline -RPROJ_ELLIPSOID.nc mars370.in -D4 -Sq0.9975 -Gmars2.nc
REM Scale to km and remove PROJ_ELLIPSOID
gmt grdmath mars.nc 1000 DIV PROJ_ELLIPSOID.nc SUB = mars.nc
gmt grdmath mars2.nc 1000 DIV PROJ_ELLIPSOID.nc SUB = mars2.nc
gmt makecpt -Crainbow -T-7/15/22 -Z > mars.cpt
gmt grdgradient mars2.nc -fg -Ne0.75 -A45 -Gmars2_i.nc
gmt grdimage mars2.nc -Imars2_i.nc -Cmars.cpt -B30g30 -BWsne -JH0/7i -P -K -X0.75i -E200 --FONT_ANNOT_PRIMARY=12p > %ps%
gmt grdcontour mars2.nc -J -O -K -C1 -A5 -Glz+/z- >> %ps%
gmt psxy -Rg -J -O -K -Sc0.045i -Gblack mars370.in  >> %ps%
echo 0 90 b) | gmt pstext -R -J -O -K -N -D-3.5i/-0.2i -F+f14p,Helvetica-Bold+jLB >> %ps%
gmt grdgradient mars.nc -fg -Ne0.75 -A45 -Gmars_i.nc
gmt grdimage mars.nc -Imars_i.nc -Cmars.cpt -B30g30 -BWsne -J -O -K -E200 -Y4.2i --FONT_ANNOT_PRIMARY=12p >> %ps%
gmt grdcontour mars.nc -J -O -K -C1 -A5 -Glz+/z- >> %ps%
gmt psxy -Rg -J -O -K -Sc0.045i -Gblack mars370.in  >> %ps%
gmt psscale -Cmars.cpt -O -K -Dx3.5i/-0.15i+jTC+w6i/0.1i+h -I --FONT_ANNOT_PRIMARY=12p -Bx2f1 -By+lkm >> %ps%
echo 0 90 a) | gmt pstext -R -J -O -N -D-3.5i/-0.2i -F+f14p,Helvetica-Bold+jLB >> %ps%
REM Clean up
del *.nc
del mars.cpt
del .gmt*
