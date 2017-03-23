


PSNAME="H_vel.ps"

#gmtset ELLIPSOID Clarke-1866
#gmtset PAPER_MEDIA A4
#gmtset BASEMAP_TYPE FANCY PLOT_DEGREE_FORMAT ddd:mm:ssF GRID_CROSS_SIZE_PRIMARY 0.05i

#range="-121/-109/48/58"
range="-180/180/-90/90"

#range="-131/-99/0/90"
#range="-131/-99/48/58"

# awk '{print $1, $2, $3}' crust1.0.model | surface -Gtemp.grd -I0.05/0.05 -R$range -T0 -Lld -Lud
surface< S362ANI.vsh.1000.txt  -Gtemp.grd -I0.5/0.5 -R$range -T0 -Lld -Lud
makecpt -T34/54/0.2 -Crainbow -I -D  > my_crust1.cpt
#psbasemap -JL-111/10/48/58/5i -R$range -B2/2:."Crustal 1.0":WSne -K -V -P > $PSNAME

######grdimage temp.grd -JN0/4.5i -R  -K -V > $PSNAME


#pscoast -JL -R -Dl -N2/3 -N1 -O -W1/100 -K -O -V >> $PSNAME

#psxy AB_Geo_new.txt -JL -R -W4/20/20/20 -M"X" -K -O -V >> $PSNAME

#grdcontour temp.grd -Ch_contour.cpt -J -R -W3/180/130/130 -A2.0f9 -K -O -V >> $PSNAME

#psscale -D6.2i/1.i/1.5i/0.33i -Cmy_crust1.cpt -B5 -O -V >> $PSNAME
#evince $PSNAME &

# ps2raster $PSNAME -A -Tef -V 
#gv $PSNAME

grdimage temp.grd -JN0/4.5i -P -K -Xc -Yc -Cmy_crust1.cpt > $PSNAME

#gmt grdcontour tt.nc -J -O -K -C500 -A1000+f10p,Helvetica,white -L500 -GL0/90/203/-10,175/60/170/-30,-50/30/220/-5 -Wa0.75p,white -Wc0.25p,white >> $ps
#gmt psxy -R -J -O -K tt.pol -W0.25p,green,. >> $ps
#pscoast -R -J -O -W1p -Gsteelblue -A0/1/1 -B30g30 -B+t"CRUST 1.0" >> $PSNAME

open $PSNAME