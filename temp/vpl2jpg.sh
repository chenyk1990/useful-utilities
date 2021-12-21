#!/bin/sh

epslist=`ls Fig/*.vpl`

for i in ${epslist[@]}
do 
name=${i%.vpl}
echo $name
vpconvert $name.vpl $name.jpg
done