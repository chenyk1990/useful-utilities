#!/bin/sh

epslist=`ls *.eps`

for i in ${epslist[@]}
name=${i%.png}
do
convert $name.eps $name.png
convert $name.png $name.pdf
done


