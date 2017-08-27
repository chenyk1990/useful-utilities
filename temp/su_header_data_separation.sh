#!/bin/sh

sulist=`ls *.su`

for i in ${sulist[@]}
do
    #remove ".su"
    ii=`echo $i |rev|cut -c4-|rev`
    echo $ii
    
    #remove header and preserve binary files
    #sustrip<$i >$ii.bin
    
    #get "NO"
    no=${ii:4}
    echo $no
    
    #get header keyword value
    sugethw <shot$no.su key=sx,gx output=geom >tmp1.txt
    ntr=`more tmp1.txt | wc -l`
    echo $ntr >tmp2.txt
    cat tmp2.txt tmp1.txt >header$no.txt
    
    #remove temporary file
    rm tmp1.txt
    rm tmp2.txt
done


