#!/bin/bash

# wrapper script to generate 
# 1. reference mesh
# 2. interface mesh,mesh index

wd=$PWD

tempdir=wdir_SURFACE; rm -fr $tempdir; mkdir -p $tempdir
cd $tempdir

outdir=$wd/RequiredFiles
recpdbfile=$outdir/"rec.pdb"
intatomlist=$outdir/"intatom.list"
refxyz_file=$outdir/"mesh.ref.xyz_"

# (1) create a reformatted rec pdbfile (CYX->CYS, HIE->HIS, noH); generate asa
echo "----------------------------------"
echo "> generating reformatted pdb and asa files ..."

reformattedfile=rec.reformat.pdb
awk '{if($3 !~/^H/) print $0}' $recpdbfile | sed -e 's/CYX/CYS/' -e 's/HI./HIS/'  > $reformattedfile
#naccess $reformattedfile # output is rec.asa
naccess $reformattedfile -r $PROTLIDv1HOME/protlid_auxPrograms/naccess2.1.1 -s $PROTLIDv1HOME/protlid_auxPrograms/naccess2.1.1  #Need to provide explicit paths of the vdw.radii and standard.data files for naccess to work

cp $reformattedfile $outdir
cp rec.asa $outdir

# (2) generate reference mesh using edtsurf (data_refXYZFileForAmber); use pdbfile with H
echo "----------------------------------"
echo "> generating full reference mesh ..."
run_generateFullXYZ.sh $recpdbfile $outdir

# (3) extract interface mesh from interface amberpdb and reference mesh
echo "----------------------------------"
echo "> generating closestAtom.xyz_ ..."
run_extractInterfaceMesh_fromIntAtomList.sh $reformattedfile $intatomlist $refxyz_file $outdir


cd $wd
