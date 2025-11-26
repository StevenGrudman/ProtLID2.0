#!/usr/bin/bash

naccessLoc=$PROTLIDv1HOME/protlid_auxPrograms/naccess2.1.1
queryPDB=$1
outdir=$2
homeDirName=$3
PDB=`basename $queryPDB`
PDB="${PDB%????}"

cd ${naccessLoc}
mkdir ${homeDirName}
cd ${homeDirName}

# (1) Removes all Hydrogen atoms from PDB
# (2) Changes (CYX->CYS, HIE->HIS)
awk '{if($3 !~/^H/) print $0}' $queryPDB | sed -e 's/CYX/CYS/' -e 's/HI./HIS/'  > tmp1.pdb
#naccess $reformattedfile # output is tmp.asa
${naccessLoc}/naccess tmp1.pdb -r ${naccessLoc} -s ${naccessLoc}

more tmp1.asa | awk '{if($10 >= 5) print $0}' > tmp2.asa
cp tmp2.asa $PDB.SA
mv $PDB.SA $outdir
cd ..
rm -r ${homeDirName}
