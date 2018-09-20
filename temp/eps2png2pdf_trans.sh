#!/bin/sh

epslist=`ls *.eps`

for i in ${epslist[@]}
do 
name=${i%.eps}
echo $name
convert $name.eps $name.png
convert $name.png $name.pdf
done