#!/bin/sh

epslist="a b c d e f g h i j"

for i in ${epslist[@]}
do
convert fig3$i.png fig3$i.pdf
done


