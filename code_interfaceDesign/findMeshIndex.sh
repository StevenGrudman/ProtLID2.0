# for each entry in the new xyz_ list; write out its corresponding index in the old list
# !!! assumes that the new list is a complete subset of the old list
if [ $# -ne 2 ]; then echo " newfile oldfile" ; exit; fi
newfile=$1
oldfile=$2

awk '{print FNR,$1,$2,$3}' $oldfile > told
n=1
while read line ; do

    x=`echo $line | awk '{print $1}' `    
    y=`echo $line | awk '{print $2}' `
    z=`echo $line | awk '{print $3}' `

    oldIndex=`awk '{if($2=="'$x'" && $3=="'$y'" && $4=="'$z'") print $1}' told`
    if [[ ! $oldIndex =~ [0-9] ]]; then echo "Warning!!!! Cannot find oldindex for line $n: $line";
    else
	echo $oldIndex
    fi
    n=$[n+1]
done < $newfile

rm -f told
