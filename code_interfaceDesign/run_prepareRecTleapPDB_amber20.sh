# generate amber-receptor files from apo/holo rcsb_pdb

wd=$PWD

if [ $# -ne 4 ] ; then echo "usage: [rcsb pdbfile] [rec-chain] [start_residno][end_residno]"; exit; fi
rcsbpdbfile=$1; 
echo $rcsbpdbfile
#if [[ $rcsbpdbfile =~ "^." ]]; then rcsbpdbfile=$wd/$rcsbpdbfile; fi
rcsbpdbfile=$wd/$rcsbpdbfile;
chainR=$2 
startR=$3
endR=$4

outdir=$wd/RequiredFiles; mkdir -p $outdir

echo $rcsbpdbfile
# check that input is pdbfile from rcsb
if [ -z "`grep "HEADER" $rcsbpdbfile`" ]; then echo "error: $rcsbpdbfile does not look like a rcsb pdbfile; quit"; exit; fi

# rcsb_index = amber_index + offset 
offset=$[startR-1] 

tempdir=wdir_AMBER; 
mkdir -p $tempdir; cd $tempdir

# extract receptor chain, remove H and OXT
grep ^ATOM $rcsbpdbfile | awk '{c=substr($0,22,1);r=substr($0,23,4)+0; if (c=="'$chainR'" && r>='$startR' && r<='$endR' && $3 !~/^H|OXT/) print $0}' > inpdbfile.tmp

#grep ^ATOM $wd/$rcsbpdbfile | awk '{c=substr($0,22,1);r=substr($0,23,4)+0; if (c=="'$chainR'" && r>='$startR' && r<='$endR' && $3 !~/^H|OXT/) print $0}' > inpdbfile.tmp

# make sure all input pdbs are checked for resid no continuity
error=`checkPDB.pl inpdbfile.tmp | wc -l`;
if [ $error -gt 0 ]; then echo "Error! inpdbfile.tmp has discontinuous residue no. Pls patch manually (CE) and so we can use patched verion. Skipping."; cd $wd; exit; fi

# get sequence and prepare disulphide bridges (1. convert CYS to CYX; 2. explicitly specify disulpfide bond)
awk '{if($1=="SSBOND" && $4=="'$chainR'" &&  ! $7=="'$chainR'") print $0}' $rcsbpdbfile > error.tmp
awk '{if($1=="SSBOND" && !$4=="'$chainR'" &&  $7=="'$chainR'") print $0}' $rcsbpdbfile >> error.tmp
if [ -s error.tmp ]; then 
    echo "Note: $pdb.$chainR has interchain disulphide bridge. quit."; cat error.tmp; rm -f error.tmp; 
    quit; 
fi
awk '{if($1=="SSBOND" && $4=="'$chainR'" &&  $7=="'$chainR'") print $5"-"$8}' $rcsbpdbfile > disulphide.tmp

grep ^ATOM inpdbfile.tmp | awk '{print substr($0,23,4),substr($0,18,3)}' | uniq | awk '{print $1,$2}' > seq.tmp

tDSfile=tleapline.disulphide.tmp; rm -f $tDSfile
for dsline in `cat disulphide.tmp`; do
    dsStart=`echo $dsline | awk -F'-' '{print $1}'`
    dsEnd=`echo $dsline | awk -F'-' '{print $2}'`
    if [ $dsStart -ge $startR ] && [ $dsStart -le $endR ] && [ $dsEnd -ge $startR ] && [ $dsEnd -le $endR ]; then 
	sed -i.bak 's/^'$dsStart' CYS/'$dsStart' CYX/'  seq.tmp
	sed -i.bak 's/^'$dsEnd' CYS/'$dsEnd' CYX/'  seq.tmp
	dsStart_i1=$[dsStart-offset]
	dsEnd_i1=$[dsEnd-offset]
	echo "bond REC.$dsStart_i1.SG REC.$dsEnd_i1.SG" >> $tDSfile
    else
	echo "Warning: disulphide bond ($dsStart $dsEnd) not in ig1: $recpdb $startR $endR"
    fi
done # end dsline

# check that all CYS in disulpide bridges are converted to CYX
if ls $tDSfile > /dev/null 2>&1; then 
    nCYXexpected=`cat $tDSfile | wc -l | awk '{print 2 * $1}'`
else
    nCYXexpected=0
fi
nCYXwritten=`grep "CYX" seq.tmp | wc -l`
if [ $nCYXexpected -ne $nCYXwritten ]; then echo "Error: CYX mismatch $nCYXexpected $nCYXwritten"; exit; fi 

# Generate and run tleap (Run 0: to get the index correct)
tfile=tleap.rec.in;
echo -e "source leaprc.protein.ff19SB \nrecseq = {" > $tfile
#echo -e "source leaprc.protein.ff19SB \n loadamberprep MSE.prepi \n loadamberparams frcmod.MSE \n recseq = {" > $tfile
awk '{print $2}' seq.tmp >> $tfile
echo -e "}" >> $tfile

echo "REC = loadpdbUsingSeq inpdbfile.tmp recseq" >> $tfile
echo "savePdb REC rec.pdb0" >> $tfile
echo "quit" >> $tfile

# run tleap and check results
logfile0=tleap.rec.0.log
tleap -f $tfile > $logfile0
grep "*" $logfile0 | grep -v "mismatch";


# Generate and run tleap (Run 1: to implement disulphide bond)
tfile=tleap.rec.in;
echo -e "source leaprc.protein.ff19SB \nrecseq = {" > $tfile
#echo -e "source leaprc.protein.ff19SB \n loadamberprep MSE.prepi \n loadamberparams frcmod.MSE \n recseq = {" > $tfile
awk '{print $2}' seq.tmp >> $tfile
echo -e "}" >> $tfile

echo "REC = loadpdbUsingSeq ./rec.pdb0 recseq" >> $tfile
if [ $nCYXexpected -gt 0 ]; then 
    cat $tDSfile >> $tfile
fi
echo "saveamberparm REC rec.prm rec.rst" >> $tfile
echo "savePdb REC rec.pdb" >> $tfile
echo "quit" >> $tfile

# run tleap and check results
logfile1=tleap.rec.1.log
tleap -f $tfile > $logfile1

grep "usage"  $logfile1
tail $logfile1 | grep "Writing";
if [ ! -s "rec.pdb" ]; then 
    echo "Error: no / empty rec.pdb"; 
fi
cp rec.pdb $outdir
cp rec.prm $outdir
echo "copied rec.pdb, rec.prm to $outdir."
echo "!!!! Always check tleap.rec.0.log and tleap.rec.1.log in wdir_AMBER to ensure tleap build is correct. !!!"

cd $wd

