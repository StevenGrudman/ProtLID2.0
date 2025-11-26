#!/bin/bash

# for each hydrophobic residue, use its fa as the rec atom (if it is in on the interface)

if [ $# -ne 1 ] ;then echo "hydinterface_pdbfile"; exit; fi
hydinterface_pdbfile=$1
wd=$PWD
masterjoblist=$PROTLIDv1HOME/sharedFiles/joblist.aafa.master0; 
hydjoblist=joblist.aafa.hyd0; grep "HYD" $masterjoblist > $hydjoblist
outfile="recFAs.hyd.list"; rm -f $outfile

# get a list of residue-residno on the hyd-interface
awk '{print $4,$5}' $hydinterface_pdbfile | sort | uniq > interface.aa-no.tmp

while read line; do

    aa=`echo $line | awk '{print $1}'`
    residno=`echo $line | awk '{print $2}'`

    # do not consider cysteines on the interface
    if [[ $aa =~ "CYS" ]]; then 
	echo " *** skipping $aa $residno ***"
	continue; 
    fi

    # 1. check if residue is a hydrophobic residue
    awk '{if($1=="'$aa'") print $3}' $hydjoblist > fas.tmp
    nFas=`cat fas.tmp | wc -l`
    #if [ $nFas -eq 0 ]; then echo "skipping $aa $residno ..."; continue; fi
    if [ $nFas -gt 1 ]; then echo "error: $aa : nFas != 1"; exit; fi

    # check if the fa is found on the interface
    for fa in `cat fas.tmp | tr "," "\n"`; do 
	if [ $fa == "RC" ]; then 
	    nFound=`awk '{if($5=='$residno' && ( ( $4~/PHE|TYR/ && $3~/C[DEGZ]/) || ($4~/TRP/ && $3~/C[EZ]|CD2|CH2/)  || ($4~/PRO/ && $3~/C[ABGD]|N/) || ($4~/HI/ && $3~/CG|[CN]D|[CN]E/) ) ) print 1}' $hydinterface_pdbfile | wc -l`
	else
	    nFound=`awk '{if( $5=='$residno' && $3=="'$fa'") print 1}'  $hydinterface_pdbfile | wc -l`; 
	fi

	if [ $nFound -gt 0 ]; then  
	    echo "$residno $aa $fa" >> $outfile
	else
	    echo "note: $aa:$residno:$fa not on interface (even though $aa:$residno has other hyd interface atoms)"
	fi
    done

done < interface.aa-no.tmp

sort -o $outfile -k1n $outfile

rm -f fas.tmp interface.aa-no.tmp