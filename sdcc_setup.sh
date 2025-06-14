## get software stack from LCG in CVMFS
## this setups both ROOT and Geant4
## docs: https://lcginfo.cern.ch/
echo "Setting up LGC_105 software stack..."
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-centos7-gcc11-opt/setup.sh

## HEP_HPC: this is not available in LGC: it's a Fermilab/NOvA-specific package
## DUNE cvmfs: compiled with an older HDF5 than LGC...
## Current solution: build HEP_HPC locally on SDCC with HDF5 from LGC.
## Use env variable to point to it later in CMakeLists.txt
export hep_hpc_path="/lbne/u/mvicenzi/hep_hpc/bin"

echo "Setup completed!"
