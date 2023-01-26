#!/bin/bash
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


module load omfit/unstable
module load intel 

cd /fusion/projects/codes/pytomo/pyTomo
python pytomo.py  $@ -t 23  --lambda_solver 0.2 -u 5 -d 5 --boundary 0 -S 7 -n --positive_constrain -G    -B 100 
