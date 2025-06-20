#!/bin/bash -l
module purge
# module load defaults
export MKL_NUM_THREADS=1
export openblas_get_num_threads=1
export OMP_NUM_THREADS=1
export NUM_THREADS=1
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/fusion/projects/codes/pytomo/SuiteSparse/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/fusion/projects/codes/pytomo/SuiteSparse/CHOLMOD/Include/
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/fusion/projects/codes/pytomo/SuiteSparse/include/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/fusion/projects/codes/pytomo/SuiteSparse/lib/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/fusion/projects/codes/pytomo/SuiteSparse/CHOLMOD/Lib/
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/fusion/projects/codes/pytomo/SuiteSparse/CHOLMOD/Lib/
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/fusion/usc/opt/intel2018/compilers_and_libraries_2018.2.199/linux/compiler/lib/intel64_lin

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/fusion/usc/opt/intel2020/intelpython3/lib/

module load omfit
module load intel

echo 'Start pyTOMO'

cd /fusion/projects/codes/pytomo/
python3  pytomo.py $@
