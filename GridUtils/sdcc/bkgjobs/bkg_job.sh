#!/bin/bash

###-----------------------------------------
### DEFINITION OF PRODUCTION PARAMETERS
###----------------------------------------

# Production name for output directories
# This is used to place output logs and files
export prodname="test_campaign"

# Define how many jobs, how many files
# Jobs will be placed in the same cluster
export n_jobs=10
export n_evt_per_job=200

# Paths to FPFsim build directory
export builddir="/direct/lbne+u/${USER}/FPFSim/build"
export fpfsim="${builddir}/FPFSim"
export setup="${builddir}/../sdcc_setup.sh"
export libdict="${builddir}/libFPFClasses_rdict.pcm"

# Generator options
export twindow="200 us"
export inputbkg="${builddir}/../backgrounds/background_input.root"

# Path to geometry macro
export geometry="${builddir}/macros/geometry_options/FLArE_only_BabyMIND.mac"
export rock="false"

# Path to log/output directory
export outdir="/gpfs01/lbne/users/fpf/$USER/CONDOR_OUTPUT"
export logdir="/gpfs01/lbne/users/fpf/$USER/CONDOR_LOGS"

###------------------------------------------------------------------
###------------------------------------------------------------------

function generate_macros {
  
  list=$1

  # for each job
  for i in $(seq 0 $((${n_jobs}-1)))
  do

    # setup different seeds per jobs
    seed1=$((42+i))
    seed2=$((47+i))

    # select the event numbering offset for this job
    istart=$((n_evt_per_job*i))

    # select name for output file 
    outputfile="${prodname}_${i}.root"
    
    # final macro path  
    macfile="${logdir}/${prodname}/mac/${prodname}_${i}.mac"
    if [ -f ${macfile} ]; then
      echo "${macfile} exists, delete file first"
      rm ${macfile}
    fi

    # add macro path to list
    echo $macfile >> $list
    rm temp.mac
    cat << EOF >> temp.mac
/control/execute ${geometry}
/det/enableRock ${rock}
/det/flare/useNativeG4Scorer true

/random/setSeeds ${seed1} ${seed2}
/run/initialize

/gen/select background
/gen/bkg/backgroundInput ${inputbkg}
/gen/bkg/backgroundWindow ${twindow}
/gen/bkg/eventOffset ${istart}

/histo/save3DEvd false
/histo/save2DEvd false
/histo/saveHit false
/histo/addDiffusion false
/histo/fileName ${outputfile}

/run/beamOn ${n_evt_per_job}
EOF
    
    # move the macro to final destination
    cp temp.mac ${macfile}
  done

}

###------------------------------------------------------------------
###------------------------------------------------------------------

function generate_submission_file {

  sub=$1; shift
  listpath=$1; shift

  cat << EOF >> ${sub}
universe                = vanilla
notification            = never
executable              = batch_bkg_script.sh
arguments               = \$(ClusterId) \$(ProcId) ${fpfsim} ${listpath} ${setup}
initialdir              = ${outdir}/${prodname}
output                  = ${logdir}/${prodname}/out/\$(ClusterId).\$(ProcId).out
error                   = ${logdir}/${prodname}/err/\$(ClusterId).\$(ProcId).err
log                     = ${logdir}/${prodname}/log/\$(ClusterId).\$(ProcId).log
getenv                  = True
request_memory          = 8000
queue ${n_jobs}
EOF

}

###------------------------------------------------------------------
###------------------------------------------------------------------

# create log directories
echo "Creating log directories for: '$prodname'"
mkdir -p ${logdir}/${prodname}/out
mkdir -p ${logdir}/${prodname}/err
mkdir -p ${logdir}/${prodname}/log
mkdir -p ${logdir}/${prodname}/mac

# check if production with the same name already exists
if [ -d "${outdir}/${prodname}" ] && [ "$(ls -A ${outdir}/${prodname})" ]; then
  echo "WARNING: ${outdir}/${prodname} is not empty!"
  echo "Outputs for production ${prodname} already exist, delete them first"
  echo "or change the name for this production!"
  exit 1
else
  mkdir -p ${outdir}/${prodname}
  if [ -d "${logdir}/${prodname}/out/" ] && [ "$(ls -A ${logdir}/${prodname}/out/)" ]; 
    then rm ${logdir}/${prodname}/out/*; fi
  if [ -d "${logdir}/${prodname}/err/" ] && [ "$(ls -A ${logdir}/${prodname}/err/)" ]; 
    then rm ${logdir}/${prodname}/err/*; fi
  if [ -d "${logdir}/${prodname}/log/" ] && [ "$(ls -A ${logdir}/${prodname}/log/)" ]; 
    then rm ${logdir}/${prodname}/log/*; fi
  if [ -d "${logdir}/${prodname}/mac/" ] && [ "$(ls -A ${logdir}/${prodname}/mac/)" ]; 
    then rm ${logdir}/${prodname}/mac/*; fi
fi


# generate all macros for each job
# as well as a list of their paths
echo "Generating list of .mac files..."
  
macrolist="macrolist.txt"
if [ -f ${macrolist} ]; then
  rm ${macrolist}
fi
touch $macrolist

generate_macros $macrolist

macropath="${logdir}/${prodname}/${macrolist}"
cp $macrolist ${macropath}

echo "Generating .sub file..."

subfile=${prodname}.sub
if [ -f ${subfile} ]; then
  rm ${subfile}
fi
touch ${subfile}

generate_submission_file $subfile $macropath

echo "Bookeeping submission files..."

cp $subfile ${logdir}/${prodname}

condor_submit $subfile

rm ${macrolist} ${subfile} temp.mac
