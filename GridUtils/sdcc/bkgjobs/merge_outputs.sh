#!/bin/bash

export prodname=$1

export dir="/gpfs01/lbne/users/fpf/${USER}/CONDOR_OUTPUT/${prodname}/${prodname}_*.root"
export out="/gpfs01/lbne/users/fpf/${USER}/CONDOR_OUTPUT/${prodname}/${prodname}.root"
list=""

for f in $(ls ${dir});
do
	#echo $f
	list+=" ${f}"
done

#echo $list
echo "Merging into ${out}"

if test -f "$out"; then
    echo "$out exists. Removing old file"
    rm $out
fi

hadd ${out} ${list}
