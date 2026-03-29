#!/bin/bash

cluster=$1
process=$2
fpfsim=$3
macrolist=$4
setup=$5

echo "Executing JOBID ${cluster}.${process}"

# source the environment
source $setup

# select the macro file
let num=$((${process}+1))
echo "Selecting macro from line ${num} in list:"

macropath=$(tail -n+${num} ${macrolist} | head -n1)
macro=`basename "$macropath"`
echo "$macro"

# running !!
echo "Running ${fpfsim} ${macropath}"
${fpfsim} ${macropath}

echo "Completed JOBID ${cluster}.${process}"
