#!/bin/bash -f

#setenv osy `fs sysname | awk -F\' '{print $2}'`

source /etc/profile.d/modules.sh

module purge 
module load mkl intel
 
export LD_LIBRARY_PATH=${MKL_HOME}/lib/intel64_lin

#load mencoder nd ffmpeg
export PATH=${PATH}/afs/@cell/common/soft/visualization/mencoder/ffmpeg/2.7/amd64_sles11/bin/
export PATH=${PATH}/afs/@cell/common/soft/visualization/mencoder/svn-2012-11-15/amd64_sles11/bin/
module load gnuplot
 

export LD_LIBRARY_PATH=${MKL_HOME}/lib/intel64_lin:/afs/ipp-garching.mpg.de/home/t/todstrci/SuiteSparse/lib/
export C_INCLUDE_PATH=${LD_LIBRARY_PATH}:/afs/ipp-garching.mpg.de/home/t/todstrci/SuiteSparse/lib/

export MKL_NUM_THREADS=1
module load anaconda/3/2020.02


rootdir=`dirname $0`                       # may be relative path
export PYTOMO=`cd $rootdir && pwd`  # ensure absolute path
cd $PYTOMO

python ./pytomo.py  $@

