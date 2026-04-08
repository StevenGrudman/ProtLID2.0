#!/bin/bash

#SBATCH --job-name=mdALA1
#SBATCH --partition=12cores,36cores
#SBATCH --array=1-21
#SBATCH --output=LOGFILES/slurm-%A_%a.out
#SBATCH --error=LOGFILES/slurm-%A_%a.err
#SBATCH --chdir=.
#SBATCH --export=ALL

# Run AMBER (minimization, md, minimization, pdb)

### here you list your tasks like (task1 task2 ...) separated by blanks
tasks=(rec.ALA.2012.1.rst rec.ALA.2207.1.rst rec.ALA.2310.1.rst rec.ALA.2312.1.rst rec.ALA.2415.1.rst rec.ALA.2506.1.rst rec.ALA.2936.1.rst rec.ALA.2942.1.rst rec.ALA.2947.1.rst rec.ALA.3134.1.rst rec.ALA.3349.1.rst rec.ALA.3548.1.rst rec.ALA.3557.1.rst rec.ALA.3568.1.rst rec.ALA.3643.1.rst rec.ALA.3661.1.rst rec.ALA.3771.1.rst rec.ALA.4061.1.rst rec.ALA.4063.1.rst rec.ALA.4173.1.rst rec.ALA.4248.1.rst)

starttime=$(date +"%T")

## default
jobindex=$((SLURM_ARRAY_TASK_ID - 1))
input=${tasks[$jobindex]}

## this is default, do not touch
wd=$PWD

## Load AMBER
module load amber/20

# MD parameters
### ROS: changed from 200 to 150
### nSnapshot=200
### SG. Change nSnapshot from 25 to 125 and snapshot_interval from 5 to 1
nSnapshot=125
n_intervalStep=100
snapshot_interval=1
cutoff=18 #999

nstlim=`echo $nSnapshot $n_intervalStep | awk '{print $1*$2}'`

# prepare the file names
shareddatadir=$PROTLIDv2HOME/sharedFiles
requiredfilesdir=$wd/RequiredFiles

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROBE TYPE
ligtype=UNCAPPED
base=`echo $input | sed 's/.rst//'`
aaTag=`echo $base | awk -F'.' '{print $2}'`
traj=`echo $base | awk -F'.' '{print $3}'`
runid=`echo $base | awk -F'.' '{print $4}'`
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# probe-related parameters
probeLength=`echo $aaTag | awk -F'-' '{print NF}'`
bProbeEven=`echo $probeLength | awk '{if($1%2==0) print 1; else print 0}'`

# handle different capping
if [[ $ligtype =~ "UNCAPPED" ]]; then cap_offset=0;
else
    if [[ $ligtype =~ "CAPPED" ]]; then cap_offset=1;
    else echo "Error: cannot handle ligtype : $ligtype"; exit; fi
fi

lastrec_residnum=`grep ^ATOM $wd/RequiredFiles/rec.pdb | tail -1 | awk '{print $5}'`
ligresidnum_start=$[lastrec_residnum + cap_offset + 1]
ligresidnum_end=$[ligresidnum_start + probeLength - 1]

# for probe restraint during MD
midIndex=`echo $probeLength | awk '{printf "%d\n",($1+1)/2}'` 
ligresidnum_mid=$[ligresidnum_start + midIndex - 1] 
if [ $bProbeEven -eq 1 ]; then 
    funcatom="C" # C of middle probe residue
else 
    funcatomfile=$shareddatadir/funcAtom_buildphase.table
    aaMid=`echo $aaTag | awk -F'-' '{print $'$midIndex'}'`
    funcatom=`grep "$aaMid" $funcatomfile | awk '{print $2}'` # FA of 1st probe residue
fi

# prepare the file names

scanprmdir=$wd/COM_PRM
initpdbdir=$wd/PDB_INIT
startrstdir=$wd/RST_INIT/$aaTag/runid_$runid
startrstfile=$input
comprmfile=rec.$aaTag.prm
xyz_file=mesh.ref.xyz_

restraintfile=restraints.in
outrst5dir="$wd/RST_05/$aaTag/runid_$runid"
outinfodir="$wd/INFO_SANDER/$aaTag/runid_$runid"
outtrj6dir="$wd/TRJ_06/$aaTag/runid_$runid"
outrst6dir="$wd/RST_06/$aaTag/runid_$runid"
outrst7dir_snap="$wd/RST_07_snapshots/$aaTag/runid_$runid"
outpdbdir="$wd/PDB_07/$aaTag" # put all runids in the same folder
outptrajlogdir="$wd/INFO_PTRAJ/$aaTag/runid_$runid"; 
mkdir -p $outtrj6dir $outrst6dir $outrst7dir_snap  $outrst5dir $outinfodir $outptrajlogdir $outpdbdir

# copy required files to remote directory
tempdir=/tmp/$input"_"$JOB_ID
rm -fr $tempdir; mkdir -p $tempdir
cp $shareddatadir/0*mi.in $tempdir
cp $shareddatadir/06md.in $tempdir
cp $shareddatadir/07mi.in $tempdir
cp $shareddatadir/$restraintfile $tempdir
cp $requiredfilesdir/$xyz_file $tempdir
cp $scanprmdir/$comprmfile $tempdir
cp $startrstdir/$startrstfile $tempdir

cd $tempdir
echo "$base"

# ============================    
# change the md input files accordingly
for inputfile in `ls 0[1-7]m*in`; do 
    sed -i 's/cut = [0-9.]*,/cut='$cutoff',/' $inputfile;
done


for i in 1 2 3 4; do 
    sed -i 's/X_LASTREC_RESIDNO_X/'$lastrec_residnum'/' 0"$i"mi.in;
    sed -i 's/X_LIG_RESIDNO_X@CG/'$ligresidnum_mid'@'$funcatom'/' 0"$i"mi.in;
#    sed -i 's/maxcyc = 1000/maxcyc = 10/' 0"$i"mi.in; # FOR QUICK TESTING ONLY; COMMENT OUT FOR PRODUCTION RUNS

done

sed -i 's/X_LASTREC_RESIDNO_X/'$lastrec_residnum'/' 06md.in 
sed -i 's/nstlim.*$/nstlim = '$nstlim',/' 06md.in
sed -i 's/ntwx.*$/ntwx = '$n_intervalStep'/' 06md.in # write out in intervals matching desired intervening steps
cp 0*in $outinfodir;

#==================================
# 1. MINIMIZE INITIAL STRUCTURE
#==================================
if [ -f $outrst5dir/$base.05mi.rst ]; then
    echo "using existing minimized rst file ..."
    cp $outrst5dir/$base.05mi.rst .
else
    now=$(date +"%T")
    echo "staring minimization $now"
    sander -O -i 01mi.in -o $base.01mi.out -p $comprmfile -c $startrstfile -ref $startrstfile  -x 01mi.trj -inf 01mi.info -r 01mi.rst
    sander -O -i 02mi.in -o $base.02mi.out -p $comprmfile  -c 01mi.rst -ref 01mi.rst -x 02mi.trj -inf 02mi.info -r 02mi.rst
    sander -O -i 03mi.in -o $base.03mi.out -p $comprmfile  -c 02mi.rst -ref 02mi.rst -x 03mi.trj -inf 03mi.info -r 03mi.rst
    sander -O -i 04mi.in -o $base.04mi.out -p $comprmfile  -c 03mi.rst -ref 03mi.rst -x 04mi.trj -inf 04mi.info -r 04mi.rst
    sander -O -i 05mi.in -o $base.05mi.out -p $comprmfile  -c 04mi.rst -ref 04mi.rst -x 05mi.trj -inf $base.05mi.info -r $base.05mi.rst
    
# if minimization fails, cp back the logfiles for debugging
    bMinimizeFail=0;
    if [ ! -f $base.05mi.rst ]; then bMinimizeFail=1; 
    else isNan=`grep NaN $base.05mi.rst`; 
	if [ -n "$isNan" ]; then  bMinimizeFail=2; fi
    fi
    if [ $bMinimizeFail -gt 0 ]; then 
	echo "$base Minimization failed: mode $bMinimizeFail (1=no 05mi; 2=nan) ";
	cp $base.*mi.out $outinfodir; rm -fr $tempdir; 
	exit; 
    fi

    cp $base.05mi.rst $outrst5dir
    cp $base.05mi.info $outinfodir

    now=$(date +"%T")
    echo "after 05mi: $now"

fi
#==================================
# 2. MD RUN
#==================================
# modify restraint file according to ligand-functional atom and receptor-atom positions
pdbreffile=$initpdbdir/rec.$aaTag.sample.pdb
ligFApos=`awk '{if($3=="'$funcatom'" && $5=='$ligresidnum_mid') print $2}' $pdbreffile`
recatompos=`awk '{if(FNR=='$traj') print $6}' $xyz_file`
if [[ ! $ligFApos =~ [0-9] ]] || [[ ! $recatompos =~ [0-9] ]]; then echo "error: restraints not correct: $ligFApos, $recatompos"; rm -fr $tempdir; exit; fi
sed -i 's/iat.*$/iat='$ligFApos','$recatompos',/' $restraintfile

cp $restraintfile $outinfodir

# one long MD run  
mdinputfile=06md.in
mdstartrstfile=$base.05mi.rst; 
trjfile=$base.06md.trj

trjfilesize_target0=100000 # just some big number
if ls $outtrj6dir/*trj > /dev/null 2>&1; then 
    trjfilesize_target1=`ls -lt $outtrj6dir/*trj | awk '{print $5}' | sort -nr | head -1`
else 
    trjfilesize_target1=0
fi
trjfilesize_target=`echo $trjfilesize_target0 $trjfilesize_target1 | awk '{if($1>$2) max=$1; else max=$2; print max}'`

now=$(date +"%T")
echo "starting md run $now"

sander -O -i $mdinputfile -o $base.06md.out -p $comprmfile -c $base.05mi.rst -ref $base.05mi.rst -x $trjfile -r $base.06md.rst
trjfilesize_actual_local=`ls -lt $trjfile | awk '{print $5}'`
if [ ! -f $trjfile ] || [ $trjfilesize_actual_local -lt $trjfilesize_target ]; then
    echo "No/incomplete 06md.trj : $trjfilesize_actual_local";
    cp $base.06md.out $outinfodir
    rm -rf $tempdir;
    exit;
fi

cp $trjfile $outtrj6dir
cp $base.06md.rst $outrst6dir
now=$(date +"%T")
echo "after md: $now"


#=============================================
# 3. MINIMIZE SNAPSHOTS, GENERATE PDB
#=============================================
ptrajlog=$base.ptraj.log; rm -f $ptrajlog
outpdbfile=$base."lig".pdb; rm -f $outpdbfile

for ((i=1; i<=$nSnapshot; i+=$snapshot_interval)); do

# generate rst file from trj snapshot    
    startrstfile_i=$base.06md.$i.rst
    cpptraj $comprmfile <<EOF >> $ptrajlog
    trajin $trjfile $i $i 
    trajout $startrstfile_i restart
EOF
   # ADDED BY CMA
   # Line below is not needed in AMBER 20 - CPPTRAJ
   # Output of PTRAJ appends i to the end of trajout command 
   #mv $startrstfile_i.$i $startrstfile_i
    
# minimization run
    basei="$base.07mi.$i"

    sander -O -i 07mi.in -o $basei.out -p $comprmfile  -c $startrstfile_i -ref $startrstfile_i -x $basei.trj -r $basei.rst
    if [ ! -f $basei.rst ]; then 
	echo "cannot generate $basei.rst"
	cp $basei.out $outinfodir; 
    fi
   
# convert min-rst to pdb
    cpptraj $comprmfile <<EOF >> $ptrajlog
trajin $basei.rst
trajout $base.pdb pdb
EOF

# extract only ligand coordinates from pdb
    echo $i | awk '{printf "%-10s%4d\n","MODEL",$1}' >> $outpdbfile
    # ADDED BY CMA
    #awk '{if($5>='$ligresidnum_start' && $5<='$ligresidnum_end') print $0}' $base.pdb.1 >> $outpdbfile
    # $base.pdb.1 is not needed in AMBER20 - CPPTRAJ
    awk '{if($5>='$ligresidnum_start' && $5<='$ligresidnum_end') print $0}' $base.pdb >> $outpdbfile
    echo "ENDMDL" >>  $outpdbfile

done

now=$(date +"%T")
echo "after pdb: $now"

nOutfile=`ls $base.07mi.*.rst | wc -l`
nExpected=`echo $nSnapshot $snapshot_interval | awk '{print $1/$2}'`; 
if [ $nOutfile -eq $nExpected ]; then 
    tarfile=$base.07mi.rst.tgz
else
    tarfile=$base.07mi.rst_partial.tgz
fi
tar -czf $tarfile $base.07mi.*.rst; cp $tarfile $outrst7dir_snap 
cp $outpdbfile $outpdbdir/
    
## now, remove the temp dir on the node
#cp -r $tempdir $wd
rm -rf $tempdir

endtime=$(date +"%T")
echo "$input $starttime $endtime"
