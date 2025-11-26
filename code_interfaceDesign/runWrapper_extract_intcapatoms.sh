#!/bin/bash
# given the interface atom list (residue,residno,atom), generate
# 1. interface capbable files for dist2recInt filtering
# 2. list of receptor-FA sites for receptor-based assignement during 
#    pharmacophore design

wd=$PWD
outdir=$wd/RequiredFiles

tempdir=wdir_INTERFACE; mkdir -p $tempdir; cd $tempdir

recpdbfile=$outdir/"rec.pdb"
intatomlist=$outdir/"intatom.list"


# ========================================================
# extract interface pdb based on intatomlist
# ========================================================
sed -e 's/CYX/CYS/' -e 's/HI./HIS/' $recpdbfile > pdb.reformatted.tmp
baseInterfacefile=rec.int.pdb; rm -f $baseInterfacefile
while read line; do
    residue=`echo $line | awk '{print $1}'`
    residno=`echo $line | awk '{print $2}'`
    atom=`echo $line | awk '{print $3}'`
    awk '{r=substr($0,23,4)+0; if(r=='$residno' && $4=="'$residue'" && $3=="'$atom'") print $0}' pdb.reformatted.tmp >> $baseInterfacefile
done < $intatomlist
# ========================================================
# interaction-capable pdbs for different interaction types
# ========================================================
extract_intcapatoms_HB.sh $baseInterfacefile 
extract_intcapatoms_HYD.sh $baseInterfacefile


# ========================================================
# extract receptor atoms for rec-based clustering
# ========================================================
extract_recFAs_HB.sh "rec.HB_HAD.pdb"
extract_recFAs_HYD.sh $baseInterfacefile

mv recFAs.*.list $outdir
mv rec.H*pdb $outdir

cd $wd
#rm -fr tempdir
