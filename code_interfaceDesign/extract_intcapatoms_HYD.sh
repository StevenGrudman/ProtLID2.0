#!/bin/bash

# extract interface atoms that are CAPABLE of interacting with a hydrophobic functional atom/gp
# !! use a lookup created on ushba's csu. the atomtype generated on everest's csu doesn't make sense !!

if [ $# -ne 1 ]; then echo "Usage: [full_interface_pdbfile]"; exit; fi
inputpdb=$1

# extract hydrophobic types per csu atom type assignment
sharedfiledir=$PROTLIDv2HOME/sharedFiles
allowedtypefile=$sharedfiledir/csu.hydrophobic.legit.list
atomtypelookup=$sharedfiledir/csu.atomtype.dat 

# extract pdblines for hydrophobic atom types
outpdbfile=rec.HYD.pdb; rm -f $outpdbfile

while read line; do

    atom=`echo $line | awk '{print $3}'`
    res=`echo $line | awk '{print $4}'`
    resno=`echo $line | awk '{print $5}'`;

    csuatomtype=`awk '{if($1=="'$res'" && $2=="'$atom'") print $3}' $atomtypelookup`
    if [[ $csuatomtype == "" ]]; then 
	echo "cannot find $line : to match $atom $res $resno in $atomtypelookup"; 
    fi
    bInteract=`grep "^$csuatomtype$" $allowedtypefile`
    if [ -n "$bInteract" ]; then echo "$line" >> $outpdbfile; fi

done < $inputpdb

#rm -f atomtype.tmp