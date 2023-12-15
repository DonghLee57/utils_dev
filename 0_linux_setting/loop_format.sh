#! /bin/bash

mylist=$(ls ./dir/*)
mylist=
for i in $mylist
do
 # format: ./dir/*_code
 code=$(echo $i |cut -f3 -d/ |cut -f2 -d_)
 if [ -f "$i" ]; then
  echo $code
else
 cp $i ./check.file
 fi
done            
