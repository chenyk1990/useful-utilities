#!/bin/bash
#		GMT EXAMPLE 07
#		$Id: example_07.sh 12321 2013-10-13 23:26:34Z pwessel $
#
# Purpose:	Make a basemap with earthquakes and isochrons etc
# GMT progs:	pscoast, pstext, psxy
# Unix progs:	echo, rm
#
ps=example_07.ps
gmt pscoast -R-50/0/-10/20 -JM9i -K -Slightblue -GP300/26:FtanBdarkbrown -Dl -Wthinnest \
	-B10 --FORMAT_GEO_MAP=dddF > $ps
gmt psxy -R -J -O -K fz.xy -Wthinner,- >> $ps
gmt psxy quakes.xym -R -J -O -K -h1 -Sci -i0,1,2s0.01 -Gred -Wthinnest >> $ps
gmt psxy -R -J -O -K isochron.xy -Wthin,blue >> $ps
gmt psxy -R -J -O -K ridge.xy -Wthicker,orange >> $ps
gmt psxy -R -J -O -K -Gwhite -Wthick -A >> $ps << END
-14.5	15.2
 -2	15.2
 -2	17.8
-14.5	17.8
END
gmt psxy -R -J -O -K -Gwhite -Wthinner -A >> $ps << END
-14.35	15.35
 -2.15	15.35
 -2.15	17.65
-14.35	17.65
END
echo "-13.5 16.5" | gmt psxy -R -J -O -K -Sc0.08i -Gred -Wthinner >> $ps
echo "-12.5 16.5 ISC Earthquakes" | gmt pstext -R -J -F+f18p,Times-Italic+jLM -O -K >> $ps
gmt pstext -R -J -O -F+f30,Helvetica-Bold,white=thin >> $ps << END
-43 -5 SOUTH
-43 -8 AMERICA
 -7 11 AFRICA
END
