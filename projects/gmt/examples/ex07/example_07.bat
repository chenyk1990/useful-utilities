REM		GMT EXAMPLE 07
REM
REM		$Id: example_07.bat 12820 2014-01-31 21:01:19Z remko $
REM
REM Purpose:	Make a basemap with earthquakes and isochrons etc
REM GMT progs:	pscoast, pstext, psxy
REM DOS calls:	del, echo
REM
echo GMT EXAMPLE 07
set ps=example_07.ps
gmt pscoast -R-50/0/-10/20 -JM9i -K -Slightblue -GP300/26:FtanBdarkbrown -Dl -Wthinnest -B10 --FORMAT_GEO_MAP=dddF > %ps%
gmt psxy -R -J -O -K fz.xy -Wthinner,- >> %ps%
gmt psxy quakes.xym -R -J -O -K -h1 -Sci -i,1,2s0.01 -Gred -Wthinnest >> %ps%
gmt psxy -R -J -O -K isochron.xy -Wthin,blue >> %ps%
gmt psxy -R -J -O -K ridge.xy -Wthicker,orange >> %ps%
echo -14.5 15.2 > tmp
echo -2 15.2 >> tmp
echo -2 17.8 >> tmp
echo -14.5 17.8 >> tmp
gmt psxy -R -J -O -K -Gwhite -Wthick -A tmp >> %ps%
echo -14.35 15.35 > tmp
echo -2.15 15.35 >> tmp
echo -2.15 17.65 >> tmp
echo -14.35 17.65 >> tmp
gmt psxy -R -J -O -K -Gwhite -Wthinner -A tmp >> %ps%
echo -13.5 16.5 | gmt psxy -R -J -O -K -Sc0.08i -Gred -Wthinner >> %ps%
echo -12.5 16.5 ISC Earthquakes | gmt pstext -R -J -F+f18p,Times-Italic+jLM -O -K >> %ps%
echo -43 -5 SOUTH > tmp
echo -43 -8 AMERICA >> tmp
echo -7 11 AFRICA >> tmp
gmt pstext -R -J -O -F+f30,Helvetica-Bold,white=thin tmp >> %ps%
del .gmt*
del tmp
