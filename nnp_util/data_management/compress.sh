#! /bin/bash
PYTHON='/bin/python'

mkdir comp_outcars
joblist=$(ls ./outcars/*)
for i in $joblist
do
 # outcars/OUTCAR_XXX
 code=$(echo $i |cut -f3 -d/ |cut -f2 -d_)
 $PYTHON ./compress_outcar.py $i comp_outcars/OUTCAR_$code
done
rm -r outcars
mv comp_outcars outcars
