#!/bin/bash

source /home/sinegordon/.bashrc
export LD_LIBRARY_PATH=/usr/local/cuda-8.0/lib64/
it1=1
it2=1

for var in 0
do
#mpirun -n 2 /usr/local/miniconda2/bin/python ./run.py 0.839 $it
/usr/local/miniconda2/bin/pytnon ./run_mix.py 0.825 0.412 $it1 $it2 0
/usr/local/miniconda2/bin/python ./run_mix.py 0.833 0.417 $it1 $it2 0
/usr/local/miniconda2/bin/python ./run_mix.py 0.841 0.420 $it1 $it2 0
done
