#!/bin/sh

epslist=`ls *.eps`

for i in ${epslist[@]}
do
epstopdf $i
done


