#!/bin/bash
#		GMT EXAMPLE 05
#		$Id: example_05.sh 11905 2013-07-06 02:23:09Z pwessel $
#
# Purpose:	Generate grid and show monochrome 3-D perspective
# GMT progs:	grdgradient, grdmath, grdview, pstext
# Unix progs:	echo, rm
#
ps=example_05.ps
gmt grdmath -R-15/15/-15/15 -I0.3 X Y HYPOT DUP 2 MUL PI MUL 8 DIV COS EXCH NEG 10 DIV \
	EXP MUL = sombrero.nc
echo '-5 128 5 128' > gray.cpt
gmt grdgradient sombrero.nc -A225 -Gintensity.nc -Nt0.75
gmt grdview sombrero.nc -JX6i -JZ2i -B5 -Bz0.5 -BSEwnZ -N-1+gwhite -Qs -Iintensity.nc -X1.5i \
	-Cgray.cpt -R-15/15/-15/15/-1/1 -K -p120/30 > $ps
echo "4.1 5.5 z(r) = cos (2@~p@~r/8) @~\327@~e@+-r/10@+" | gmt pstext -R0/11/0/8.5 -Jx1i \
	-F+f50p,ZapfChancery-MediumItalic+jBC -O >> $ps
rm -f gray.cpt sombrero.nc intensity.nc
