#!/bin/bash

# script to generate full SAS mesh
recpdbfile=$1
outdir=$2

wd=$PWD

# run EDTSurf
echo "running EDTSurf ..."
scalefactor=1 # ~ reciprocal of mesh length (in A)
proberadius=1 # in Angstrom
surfacetype=2 # 2=SAS
side=2 # outside
outname="rec" # can be any name

EDTSurf -s $surfacetype -i $recpdbfile -p $proberadius -h $side -f $scalefactor -o $outname > $outname.log 2>&1
rm -f rec.asa  rec.log

# extract mesh from plyfile
plyfile=$outname.ply
xyz_file=$outname.xyz_
extractPLY.pl $plyfile > $xyz_file # (x,y,z) only

# generate reference xyz that includes nearest heavy atom
reformattedFullpdbfile="rec.reformat.pdb"
awk '{if($3 !~/^H/) print $0}' $recpdbfile | sed -e 's/CYX/CYS/' -e 's/HI./HIS/'  > $reformattedFullpdbfile

refxyz_file=mesh.ref.xyz_
addClosestAtom_toXYZ.pl $xyz_file $reformattedFullpdbfile > $refxyz_file

cp $refxyz_file $outdir
cp $reformattedFullpdbfile $outdir

