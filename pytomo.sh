#!/bin/tcsh -f


setenv osy `fs sysname | awk -F\' '{print $2}'`

source /etc/profile.d/modules.csh

module purge 
module load mkl intel
setenv LD_LIBRARY_PATH ${MKL_HOME}/lib/intel64_lin
    
#load mencoder nd ffmpeg
set path = ($path /afs/@cell/common/soft/visualization/mencoder/ffmpeg/2.7/amd64_sles11/bin/)
set path = ($path /afs/@cell/common/soft/visualization/mencoder/svn-2012-11-15/amd64_sles11/bin/)
module load gnuplot/5.2

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/afs/ipp-garching.mpg.de/home/t/todstrci/SuiteSparse/lib/
setenv C_INCLUDE_PATH ${LD_LIBRARY_PATH}:/afs/ipp-garching.mpg.de/home/t/todstrci/SuiteSparse/lib/

setenv MKL_NUM_THREADS 1
module load anaconda/3/2019.03

python   /afs/ipp/home/g/git/python/tomo/pytomo.py  $argv

