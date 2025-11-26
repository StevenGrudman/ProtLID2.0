#!/bin/bash

# set up required files for ProtLID 
# input: 
# - complex's RCSB PDB file (with SSBOND specified if any), 
# - receptor and ligand regions to consider (residue number ranges) 

export AMBERHOME=/usr/local/bio/amber20

wd=$PWD
outdir="RequiredFiles"; mkdir -p $outdir

if [ $# -ne 4 ]; then echo "usage: [compdbfile] [chainR][Rstart][Rend]"; exit; fi
compdbfile=$1
chainR=$2
startR=$3
endR=$4

# 1. generate AMBER-format receptor pdb,prm using tleap 
echo "=============================="
echo "generating AMBER-format pdb,prm ..."
echo "=============================="
run_prepareRecTleapPDB_amber20.sh $compdbfile $chainR $startR $endR

# 2. extract interface PDBs based on intatomlist
echo "=============================="
echo "generating interface PDBs ... "
echo "=============================="
runWrapper_extract_intcapatoms.sh

# 3. generate interface mesh
echo "=============================="
echo "generating interface mesh ... "
echo "=============================="
runWrapper_generateIntMesh.sh
