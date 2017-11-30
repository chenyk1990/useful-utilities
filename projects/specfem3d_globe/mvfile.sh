#!/bin/bash

eventlist=`more cmtlist`

for ievent in ${eventlist[@]}
do      
rm  -rf $ievent
done