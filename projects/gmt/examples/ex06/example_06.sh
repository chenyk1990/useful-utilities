#!/bin/bash
#		GMT EXAMPLE 06
#		$Id: example_06.sh 12017 2013-08-03 18:17:24Z jluis $
#
# Purpose:	Make standard and polar histograms
# GMT progs:	pshistogram, psrose
# Unix progs:	rm
#
ps=example_06.ps
gmt psrose fractures.d -: -A10r -S1.8in -P -Gorange -R0/1/0/360 -X2.5i -K -Bx0.2g0.2 \
	-By30g30 -B+glightblue -W1p > $ps
gmt pshistogram -Bxa2000f1000+l"Topography (m)" -Bya10f5+l"Frequency"+u" %" \
	-BWSne+t"Histograms"+glightblue v3206.t -R-6000/0/0/30 -JX4.8i/2.4i -Gorange -O \
	-Y5.0i -X-0.5i -L1p -Z1 -W250 >> $ps
