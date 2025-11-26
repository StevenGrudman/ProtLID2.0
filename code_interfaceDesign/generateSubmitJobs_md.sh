#!/bin/bash

# generate submission scripts to perform md on cluster

templatedir=$PROTLIDv1HOME/code_interfaceDesign
templatefile=$templatedir/template.1-7_pdb.sh
runname=md

scriptfilebase=`basename $templatefile .sh | sed 's/template/submit/'`

for probe in ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL; do
    for runid in `ls -d RST_INIT/$probe/runid_* | sed 's/^.*_//'`; do
	ls RST_INIT/$probe/runid_$runid/*rst | sed 's/^.*\///' > newjoblist
	
	if [ ! -s newjoblist ]; then echo "$probe-$runid: no jobs"; continue; fi
	newscriptfilebase=$scriptfilebase.$probe.$runid
	$PROTLIDv1HOME/scripts/writeSubmitscript.sh $templatefile newjoblist > $newscriptfilebase.sh;
	
	sed -i 's/#$ -N .*$/#$ -N '$runname$probe$runid'/' $newscriptfilebase.sh;
	
	echo "written to $newscriptfilebase.sh"
	
    done

done
