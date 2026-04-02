#!/bin/bash


# given a complex pdb file and rec,lig residnue number range,
# generate a list of interface atoms (residue,residno,atom) 
# using the following interface definition:
# 
# - INTERCAAT dist <= 4 
# - ASA > cutoff
# - all hbond atoms (i.e. all N,O element)
# - all S/C centers of hydrophobic residues

wd=$PWD

if [ $# -ne 7 ]; then echo "usage: [compdbfile] [chainR][Rstart] [Rend] [chainL][Lstart][Lend]"; exit; fi
compdbfile=$1; if [[ $compdbfile =~ "^." ]]; then compdbfile=$wd/$compdbfile; fi
chainR=$2
startR=$3
endR=$4
chainL=$5
startL=$6
endL=$7

asaCutoff=1 # min exposed SA (in Angstrom^2) for an atom to be considered interface

sharedfilesdir=$PROTLIDv2HOME/sharedFiles
masterjoblist=$sharedfilesdir/joblist.aafa.master0 
outdir=$wd/RequiredFiles; mkdir -p $outdir

mkdir -p tempdir; cd tempdir

# extract receptor chain, renumber residues starting from 1
grep ^ATOM ../$compdbfile | awk '{c=substr($0,22,1);r=substr($0,23,4)+0; if (c=="'$chainR'" && r>='$startR' && r<='$endR' && $3 !~/^H/) print $0}' > rec.pdb
renumberPDBResidues.pl -i rec.pdb; mv tmp_rec.pdb rec.pdb

# extract ligand chain
grep ^ATOM ../$compdbfile | awk '{c=substr($0,22,1);r=substr($0,23,4)+0; if (c=="'$chainL'" && r>='$startL' && r<='$endL' && $3 !~/^H/) print $0}' > lig.pdb

# ================================
# 1) get interaction per INTERCAAT
# ================================
echo "generating intercaat interactions ..."

dist_cutoff=4

intercaatfile=com.intercaat

cat rec.pdb > com.pdb
echo "TER"   >> com.pdb
cat lig.pdb >> com.pdb
echo "TER"   >> com.pdb
echo "END"  >> com.pdb

cd ..
python $PROTLIDv2HOME/protlid_auxPrograms/intercaat/intercaat.py -fp ./tempdir/ -qc $chainR -ic $chainL -di no -pdb com.pdb > $intercaatfile
mv com.intercaat tempdir
cd tempdir
# remove first line
sed -i '1d' com.intercaat

# filter by distance
awk -v distcutoff=$dist_cutoff '{if($11<=distcutoff) print $0}' $intercaatfile > temp.out;  

awk '{print $1,$2,$4}' temp.out | sort -k2 -n | uniq > intercaat_intatom.list # residue,residno,atom

# ============================================
# 2) get solvent exposed atoms using NACCESS
# ============================================
echo "generating solvent accessible file ..." 
asagtXfile=rec.sagt$asaCutoff.asa
#naccess rec.pdb # output is rec.asa
naccess rec.pdb -r $PROTLIDv2HOME/protlid_auxPrograms/naccess2.1.1 -s $PROTLIDv2HOME/protlid_auxPrograms/naccess2.1.1
awk '{a=substr($0,55,8)+0;if(a>'$asaCutoff') print $0}' rec.asa > $asagtXfile
# ========================================================
# 3) extract interface atoms per intercaat-list from asagtXfile, 
#    only extract S or C from hydrophobic residues
# ========================================================
outfile=$outdir/intatom.list; rm -f $outfile

grep HYD $masterjoblist | awk '{print $1}' | sort | uniq > hydres.list

echo "ready to write to $outfile"
while read line; do
    residue=`echo $line | awk '{print $1}'`
    residno=`echo $line | awk '{print $2}'`
    atom=`echo $line | awk '{print $3}'`

    # do not consider cysteines on the interface
    if [[ $residue =~ "CYS" ]]; then 
        echo " *** skipping $residue $residno ***"
        continue; 
    fi

    # check if this residue is hydrophobic
    hydFound=`grep $residue hydres.list`; if [ -n "$hydFound" ]; then bHydRes=1; else bHydRes=0; fi

    bExtract=`echo $atom $residue $bHydRes | awk '{if($1 ~/^O/ || ($1~/^N/ && $2 !="PRO") || ($1~/^S|^C/ && $3==1) ) print 1; else print 0}'`
    if [ $bExtract -eq 0 ]; then continue; fi

    echo "$residno -> $[residno-pdboffset] $residue $atom"
                                                                                 
    bResMismatch=`awk -v no=$residno -v ue=residue  '{r=substr($0,23,4)+0; if(r==no && ! $4==ue) print 1}' $asagtXfile | wc -l`
    if [ $bResMismatch -gt 0 ]; then 
        echo "the offset given causes a mismatch between pdbfile and intercaat reslist:" >&2
        echo "intercaatline: $line" >&2
        echo "pdbline: " >&2
        awk -v no=$residno -v ue=residue '{r=substr($0,23,4)+0; if(r==no && $4==ue) print $0}' $asagtXfile | head -1 >&2
        exit
    fi
      
    awk -v no=$residno -v ue=$residue -v om=$atom '{r=substr($0,23,4)+0; if(r==no && $4==ue && $3==om) print $4,r,$3}' $asagtXfile > t
    if [ ! -s t ]; then echo "warning: did not find matching residues for intercaat $line => $residno"; continue; fi
    cat t >> $outfile

done < intercaat_intatom.list

sort -o $outfile -k2n $outfile
echo "written to $outfile"
cd $wd
#rm -fr tempdir

