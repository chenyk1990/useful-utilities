#!/bin/bash
#		GMT EXAMPLE 31
#		$Id: example_31.sh 15178 2015-11-06 10:45:03Z fwobbe $
#
# Purpose:	Illustrate usage of non-default fonts in PostScript
# GMT progs:	gmtset, pscoast, psxy, pstext, pslegend
# Unix progs:	gs, awk, cat, rm
#
file=example_31
ps=${file}.ps
ps_outlined=${file}_outlined.ps
eps_outlined=${file}_outlined.eps

# create file PSL_custom_fonts.txt in current working directory
# and add PostScript font names of Linux Biolinum and Libertine
$AWK '{print $1, 0.700, 0}' << EOF > PSL_custom_fonts.txt
LinBiolinumO
LinBiolinumOI
LinBiolinumOB
LinLibertineOB
EOF

# common settings
gmt gmtset FORMAT_GEO_MAP ddd:mm:ssF \
MAP_DEGREE_SYMBOL colon \
MAP_TITLE_OFFSET 20p \
MAP_GRID_CROSS_SIZE_PRIMARY 0.4c \
PS_LINE_JOIN round \
PS_CHAR_ENCODING ISO-8859-5 \
FONT LinBiolinumO \
FONT_TITLE 24p,LinLibertineOB \
MAP_ANNOT_OBLIQUE 42

# map of countries
gmt pscoast -Dl -R-7/31/64/66/r -JL15/50/40/60/16c -P \
	-Bx10g10 -By5g5 -B+t"Europe\072 Countries and Capital Cities" -A250 \
	-Slightblue -Glightgreen -W0.25p -N1/1p,white -K > $ps

# mark capitals
gmt psxy europe-capitals-ru.csv -R -J -i0,1 \
-Sc0.15c -G196/80/80 -O -K >> $ps

# small EU cities
$AWK 'BEGIN {FS=","} $4 !="" && $4 <= 1000000 {print $1, $2}' europe-capitals-ru.csv | \
gmt psxy -R -J -Sc0.15c -W0.25p -O -K >> $ps

# big EU cities
$AWK 'BEGIN {FS=","} $4 > 1000000 {print $1, $2}' europe-capitals-ru.csv | \
gmt psxy -R -J -Sc0.15c -W1.25p -O -K >> $ps

# label big EU cities
$AWK 'BEGIN {FS=","} $4 > 1000000 {print $1, $2, $3}' europe-capitals-ru.csv | \
gmt pstext -R -J -F+f7p,LinBiolinumOI+jBL -Dj0.1c -Gwhite -C5% -Qu -TO -O -K >> $ps

# construct legend
cat << EOF > legend.txt
G -0.1c
H 10 LinBiolinumOB Population of the European Union capital cities
G 0.15c
N 2
S 0.15c c 0.15c 196/80/80 0.25p 0.5c < 1 Million inhabitants
S 0.15c c 0.15c 196/80/80 1.25p 0.5c > 1 Million inhabitants
N 1
G 0.15c
L 8 LinBiolinumOB L Population in Millions 
N 6
EOF

# append city names and population to legend
$AWK 'BEGIN {FS=","; f="L 8 LinBiolinumO L"}
  $4 > 1000000 {printf "%s %s:\n%s %.2f\n", f, $3, f, $4/1e6}' \
  europe-capitals-ru.csv >> legend.txt

# reduce annotation font size for legend
gmt gmtset FONT_ANNOT_PRIMARY 8p

# plot legend
gmt pslegend -R -J -DjTR+o0.1c+w8.0c+l1.2 \
-C0.3c/0.4c -F+p+gwhite -O legend.txt >> $ps

# make a PostScript and a PDF file with outlined fonts
# unfortunately gmt psconvert won't be able to crop that file correctly anymore
# use Heiko Oberdiek's pdfcrop (http://code.google.com/p/pdfcrop2/) instead
# or crop with gmt psconvert -A -Te before
#
# a. remove GMT logo and crop EPS:
#gmt psconvert -P -Au -Te -C-sFONTPATH="${PWD}/fonts" -Fex31CropNoLogo $ps
# b. make PS with outlined fonts:
#gs -q -sPAPERSIZE=a3 -dNOCACHE -dSAFER -dNOPAUSE -dBATCH -dNOPLATFONTS \
#  -sDEVICE=pswrite -sFONTPATH="${PWD}/fonts" -sOutputFile=$ps_outlined ex31CropNoLogo.eps
# c. make croppepd EPS:
#gs -q -dNOCACHE -dSAFER -dNOPAUSE -dBATCH -dEPSCrop -sDEVICE=epswrite \
#  -sOutputFile=$eps_outlined $ps_outlined
# d. make cropped PDF:
#gmt psconvert -P -A -Tf $ps_outlined
# uncomment to do conversation to PDF and PNG
# you will get a PDF with subsetted TrueType/PostScript fonts embedded
# which you can still edit with your favorite vector graphics editor
#export GS_FONTPATH="${PWD}/fonts"
#gmt psconvert -P -A -Tf $ps
#gmt psconvert -P -A -Tg -E110 $ps
# clean up
rm -f gmt.history gmt.conf PSL_custom_fonts.txt legend.txt ex31CropNoLogo.eps

exit 0
