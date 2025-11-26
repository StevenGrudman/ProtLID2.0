# extract interface mesh based on given interface file and full reference mesh
reformattedFullpdbfile=$1
intatomlist=$2
refxyz_file=$3
outdir=$4

wd=$PWD

# ======== generate interface mesh based on interface atom list =======================
outmeshfile=$outdir/"mesh.xyz_"
findClosestXYZs_byClosestAtomFullStruc.pl $refxyz_file $reformattedFullpdbfile $intatomlist > $outmeshfile
echo "written to $outmeshfile"

# ======== figure out the corresponding indices in the reference xyz_file =======================
n=`cat $outmeshfile | wc -l  | awk '{print $1}'`
outindexfile=$outdir/"meshindex.$n.list"
findMeshIndex.sh $outmeshfile $refxyz_file > $outindexfile
echo "written to $outindexfile"
