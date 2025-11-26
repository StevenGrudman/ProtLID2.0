#!/bin/bash

# extract list of hydrogen bond atoms (residno,aa,atom) from interface file
# for use in receptor-based assignment during pharmacophore building

if [ $# -ne 1  ] ;then echo "hbinterface_pdbfile"; exit; fi
interfacefile_hb=$1
outfile="recFAs.hbond.list"; rm -f $outfile

while read line; do
    
    aa=`echo $line | awk '{print $4}'`
    residno=`echo $line | awk '{print $5}'`
    fa=`echo $line | awk '{print $3}'`
	echo "$residno $aa $fa" >> $outfile; 

done <$interfacefile_hb
