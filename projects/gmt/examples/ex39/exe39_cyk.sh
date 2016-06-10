#!/bin/bash
#               GMT EXAMPLE 39
#               $Id$
#
# Purpose:      Illustrate evaluation of spherical harmonic coefficients
# GMT progs:    psscale, pstext, makecpt, grdimage, grdgradient, sph2grd
# Unix progs:   rm
#
ps=example_39_t10.ps

# Evaluate the first 180, 90, and 30 order/degrees of Venus spherical
# harmonics topography model, skipping the L = 0 term (radial mean).
# File truncated from http://www.ipgp.fr/~wieczor/SH/VenusTopo180.txt.zip
# Wieczorek, M. A., Gravity and topography of the terrestrial planets,
#   Treatise on Geophysics, 10, 165-205, doi:10.1016/B978-044452748-6/00156-5, 2007

sph2grd VenusTopo180.txt -I1 -Rg -Ng -Gv2.nc -F1/1/85/90
#psscale --FORMAT_FLOAT_MAP="%'g" -Ct.cpt -K -Dx1.25i/-0.2i+jTC+w5.5i/0.1i+h -Bxaf -By+lm > $ps
#grdgradient v2.nc -Nt0.75 -A45 -Gvint.nc

#grdimage v2.nc -Ivint.nc -JG -O -K -Bg -Ct.cpt -Xc -Yc >> $ps

#grdimage v2.nc -JN0/4.5i -B  -K -Xc -Yc >$ps
#grdimage v2.nc -J -B  >$ps

grdimage v2.nc -JN0/4.5i -Xc -Yc >$ps


rm -f v*.nc t.cpt

open $ps 