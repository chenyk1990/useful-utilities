#!/bin/sh

epslist=`ls *.pdf`

for i in ${epslist[@]}
do
name=${i%.pdf}
convert $name.pdf $name.eps
done


