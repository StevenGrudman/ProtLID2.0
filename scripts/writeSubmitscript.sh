#!/bin/sh

if [ $# -ne 2 ]; then echo "Usage: [script template] [joblist to insert]"; exit; fi
 
scriptfile=$1
jobfile=$2

# copy the top half of scriptfile, until tasks=()
sed -n '1,/^tasks/ p' $scriptfile | sed '/^tasks/d' > script.tmp

# put in the tasks
cat $jobfile | tr "\n" " " | sed -e 's/^/tasks=(/' -e 's/$/)\n/' >> script.tmp
ntask=`cat $jobfile | wc -l`;

# copy the rest of the scriptfile
sed -n '/^tasks/,$ p' $scriptfile | sed '/^tasks/d' >> script.tmp

# change the number of jobs
sed -i 's/^#$ -t .*$/#$ -t 1-'$ntask'/' script.tmp

cat script.tmp

rm -f script.tmp
