#!/bin/bash

# extract interface atoms that are CAPABLE of interacting with a hydrogen HD/HA or HAD
# 
if [ $# -ne 1 ]; then echo "Usage: [full interface file]"; exit; fi
inputpdbfile=$1
hbondatomdeffile=$PROTLIDv1HOME/sharedFiles/def.Hbond.dat

# list atoms that are HD-, HA-capable
for type in HA HD; do
    hb_aafa_outfile=aafa.$type.list; rm -f $hb_aafa_outfile
    for jobline in `awk '{if($4=="'$type'" || $4=="HAD") print $1":"$3}' $hbondatomdeffile`; do
	aa=`echo $jobline | awk -F':' '{print $1}'`
	for fa in `echo $jobline | awk -F':' '{print $2}' | tr "," "\n"`; do
	    echo $aa:$fa >> $hb_aafa_outfile
	done # fa
    done # jobline

done # type

# for each of the interface residue, extract pdblines from atoms capable for forming HA or HD
# HAD capable=union of HA-capable and HD-capable atoms
for type in HA HD; do
    
    outpdbfile=rec.HB_"$type".pdb; rm -f $outpdbfile
    hb_aafa_outfile=aafa.$type.list
    for jobline in `cat $hb_aafa_outfile`; do
	aa=`echo $jobline | awk -F':' '{print $1}'`
	fa=`echo $jobline | awk -F':' '{print $2}'`
	awk '{if( $4=="'$aa'" && $3=="'$fa'") print $0}' $inputpdbfile >> $outpdbfile
    done

    # include backbone N and O respectively
    awk '{if( "'$type'"=="HD" && $3=="N" && $4 !~/PRO/ ) print $0}' $inputpdbfile >> $outpdbfile
    awk '{if( "'$type'"=="HA" && $3=="O") print $0}' $inputpdbfile >> $outpdbfile

    sort -k5n $outpdbfile | uniq > t; mv t $outpdbfile
done # HA,HD

# HAD is just a union of HA and HD capable atoms
outpdbfile=rec.HB_HAD.pdb
cat rec.HB_HA.pdb rec.HB_HD.pdb | sort | uniq > $outpdbfile
sort -k5n $outpdbfile -o $outpdbfile

rm -f aafa.HA.list aafa.HD.list