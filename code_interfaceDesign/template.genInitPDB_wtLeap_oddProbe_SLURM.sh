#!/bin/bash

#SBATCH --job-name=tleap
#SBATCH --partition=36cores
#SBATCH --array=1-19
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.out
#SBATCH --chdir=.
#SBATCH --export=ALL
#SBATCH --requeue

#  ******* FOR PROBES WITH ODD NUMBERS OF RESIDUES *********
#
# generate probes initial coordinates using Amber's tleap + MODELLER
# by specifying the position of
# 1. N (amino nitrogen of middle residue)
# 2. C (carboxyl carbon of middle residue)
# 3. FA1 (side chain atom near terminal)
# 4. FA2 (atom next to FA1, proximal to main-chain)
# - rotate N-C orientation for different runs
#

# generate probes initial coordinates using Amber's tleap + MODELLER
# - orient the probe's FA inwards and mainchain outwards
# - rotate the N-C orientation for different runs

### here you list your tasks like (thisfile1.perl thisfile2.perl ... thisfile1456.perl)
### separated by blanks

tasks=(ALA ARG ASN ASP CYS GLN GLU HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL)
#tasks=(ALA-PHE-ALA)

## default
jobindex=$((SLURM_ARRAY_TASK_ID - 1))
input=${tasks[$jobindex]}


wd=$PWD

# executables
modeller=modeller # symbolic link to e.g. mod9.14,mod9.15 etc so we don't have to hardcode the version

# directories
shareddir=$PROTLIDv2HOME/sharedFiles
requiredfilesdir=$wd/RequiredFiles
interfaceDesign_scriptdir=$PROTLIDv2HOME/code_interfaceDesign
# specify files
#~~~~~~~~~~~
aaTag=$input
nRun=7
#~~~~~~~~~~~
aaSeq_threeLetter=`echo $aaTag | sed 's/-/ /g'` # e.g. ALA PHE ALA
probeLength=`echo $aaTag | awk -F'-' '{print NF}'`
midIndex=`echo $probeLength | awk '{print ($1+1)/2}'` # for odd-probe: middle 1 residue
aaMid=`echo $aaTag | awk -F'-' '{print $'$midIndex'}'` 

funcfile=$shareddir/funcAtom_buildphase.table

tfile=tleap.$aaTag.com.in
recpdbfile=rec.pdb
meshindexfile=`ls $requiredfilesdir/meshindex.[1-9]*.list | sed 's/^.*\///'`
xyz_file=mesh.ref.xyz_

amberdir=$wd/wdir_AMBER
scantleapdir=$wd

outprmdir=$wd/COM_PRM
outrstdir=$wd/RST_INIT/$aaTag
outpdbdir=$wd/PDB_INIT
mkdir -p $outprmdir $outrstdir $outpdbdir

# copy files to remote dir
tempdir=/tmp/"$input"_$JOB_ID
mkdir -p $tempdir
cp $shareddir/initial_probe_plus_receptor_build.py $tempdir
cp $shareddir/probe_plus_receptor_build.py $tempdir
cp $amberdir/tleap.rec.in $tempdir
cp $requiredfilesdir/$recpdbfile $tempdir
cp $requiredfilesdir/rec.reformat.pdb $tempdir
cp $requiredfilesdir/$xyz_file $tempdir
cp $requiredfilesdir/$meshindexfile $tempdir
cp $interfaceDesign_scriptdir/computePoints_oddProbe.pl $tempdir
cp $interfaceDesign_scriptdir/generatePIR.pl $tempdir

cd $tempdir

#=====
# determine the last residue no of the receptor  
lastRecResidno=`grep ^ATOM rec.reformat.pdb | tail -1 | awk '{print $5}'`

# side-chain atoms to be restrained
fa1=`grep "$aaMid" $funcfile | awk '{print $2}'`
fa2=`grep "$aaMid" $funcfile | awk '{print $3}'`
if [ -z "$fa1" ] || [ -z "$fa2" ]; then echo "error: cannot determine fa1 ($fa1) or fa2 ($fa2). quit"; rm -fr $tempdir; exit; fi
echo "fa1 = $fa1"
echo "fa2 = $fa2"

# generate tleap.com.in file 
grep -v "^save\|quit" tleap.rec.in | sed 's/rec.pdb0/rec.pdb/' > $tfile 
echo -e "
ligSeq = {$aaSeq_threeLetter}
LIG = loadpdbUsingSeq ./temp.pdb ligSeq
saveamberparm LIG lig.$aaTag.prm lig.$aaTag.rst
savePdb LIG lig.$aaTag.pdb
COM = combine {REC LIG}
saveamberparm COM rec.$aaTag.prm rec.$aaTag.0.rst
savePdb COM rec.$aaTag.0.pdb
quit
" >> $tfile

# generate ali file for Modeller
./generatePIR.pl $aaTag > temp_alignment.ali 

# ====================================================
# generate multiple probe-position at each mesh point
# ====================================================
for n in `cat $meshindexfile`; do

    line=`sed -n ''$n'p' $xyz_file`

# X Y Z position
    xm=`echo $line | awk '{print $1}'`
    ym=`echo $line | awk '{print $2}'`
    zm=`echo $line | awk '{print $3}'`
    recatom_index=`echo $line | awk '{print $6}'`

# given mesh and coords of the closest rec atom, then compute PN,P2,PC in multiple orientations
    bondlength=1.55
    xyzrecatom=`awk '{if($2=='$recatom_index') print $6,$7,$8}' $recpdbfile`
    unitVecR_FA=`echo $xm $ym $zm "$xyzrecatom" | awk '{dx=$1-$4;dy=$2-$5;dz=$3-$6;n=sqrt(dx^2+dy^2+dz^2); print dx/n,dy/n,dz/n}'`
    ./computePoints_oddProbe.pl $xm $ym $zm $xyzrecatom $bondlength $nRun > positions.tmp
    r=1 # runid
    while read line; do

	angle=`echo $line | awk '{print $1}'`
	xN=`echo $line | awk '{print $2}'`
	yN=`echo $line | awk '{print $3}'`
	zN=`echo $line | awk '{print $4}'`
	x2=`echo $line | awk '{print $5}'`
	y2=`echo $line | awk '{print $6}'`
	z2=`echo $line | awk '{print $7}'`
	xC=`echo $line | awk '{print $8}'`
	yC=`echo $line | awk '{print $9}'`
	zC=`echo $line | awk '{print $10}'`
	
	# generate ligand pdbfile on the fly 
	echo "N    $aaMid $xN $yN $zN $midIndex" | awk '{printf "%-6s %4d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f %5.2f %5.2f\n","ATOM",1,$1,$2," ",$6,$3,$4,$5,1.0,0.0}' >  temp_centerRes_reduced.pdb
	echo "$fa1 $aaMid $xm $ym $zm  $midIndex" | awk '{printf "%-6s %4d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f %5.2f %5.2f\n","ATOM",1,$1,$2," ",$6,$3,$4,$5,1.0,0.0}' >> temp_centerRes_reduced.pdb
	echo "$fa2 $aaMid $x2 $y2 $z2 $midIndex" | awk '{printf "%-6s %4d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f %5.2f %5.2f\n","ATOM",1,$1,$2," ",$6,$3,$4,$5,1.0,0.0}' >> temp_centerRes_reduced.pdb
	echo "C    $aaMid $xC $yC $zC $midIndex" | awk '{printf "%-6s %4d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f %5.2f %5.2f\n","ATOM",1,$1,$2," ",$6,$3,$4,$5,1.0,0.0}' >> temp_centerRes_reduced.pdb

###################################################################################################################################################
	# USE MODELLER TO BUILD PROBE (based on Nelson's code)
#####################################################################################################################################################

	# run Modeller first without minimization and default restraints to get .ini file (to know atom indices in order to generate restraint file)
	$modeller initial_probe_plus_receptor_build.py # this script exits after producing initial model

	echo "==============="

	# generate reduced probe restraint file
	grep ^ATOM PROBE_PLUS_RECEPTOR.ini | awk '{ if ($6=='$[lastRecResidno+midIndex]' && ($3=="'$fa1'" || $3=="'$fa2'")) {for(i=0;i<3;i++) { print "R", 3, 1, 9+i, 34, 1, 2, 0, $2, $(7+i), 0.0001;} } }' > REDUCED_PROBE_RESTRAINT.rsr
	
	# generate receptor restraint file
	grep ^ATOM PROBE_PLUS_RECEPTOR.ini | awk '{ if ($6<='$lastRecResidno') {for(i=0;i<3;i++) { print "R", 3, 1, 9+i, 34, 1, 2, 0, $2, $(7+i), 0.0001;} } }' > RECEPTOR_POSITION_RESTRAINT.rsr

	# run Modeller appending the restraints just generated to the default ones to build probe and receptor
	$modeller probe_plus_receptor_build.py 
	
	# extract probe from Modeller model
	grep ^ATOM PROBE_PLUS_RECEPTOR.B99990001.pdb |  awk '{if($6>'$lastRecResidno' && $3!="OXT") print $0}' > temp.pdb

#####################################################################################################################################################
	# END MODELLER
#####################################################################################################################################################

	# specify index for rst file 
	sed -i.bak  's/\.0.rst/.'$n.$r'.rst/' $tfile
	sed -i.bak  's/\.[0-9.-]*.rst/.'$n.$r'.rst/' $tfile

	sed -i.bak  's/\.0.pdb/.'$n.$r'.pdb/' $tfile
	sed -i.bak  's/\.[0-9.-]*.pdb/.'$n.$r'.pdb/' $tfile
	
	# run tleap
	tleap -s -f $tfile > log2 #/dev/null 2>&1
	r=$[r+1]

    done < positions.tmp
done 

cat log log2 > $aaTag.log
#=================
# copy back files
#chmod -R 755 $tempdir

samplecompdb=`ls rec.$aaTag.*.pdb | head -1`; cp $samplecompdb $outpdbdir/rec.$aaTag.sample.pdb

for r in `seq 1 $nRun`; do 
    mkdir -p $outrstdir/runid_$r; 
    cp rec.$aaTag.*.$r.rst $outrstdir/runid_$r;
done

cp $tfile $outprmdir
cp rec.*.prm $outprmdir
cp lig.*.prm $outprmdir
cp $aaTag.log $outprmdir


# remove tmp dir
#cp -r $tempdir $wd
rm -fr $tempdir

#cd $wd
