#!/bin/bash
# Название расчитываемой задачи. Может быть любым.
#SBATCH --job-name="DMPC"
#
# Множество вычислительных узлов для расчета задачи. Определяет характеристику
# вычислительных узлов.
#SBATCH --partition=intelv4-batch
#
# Запускать каждый расчет на одном узле.
#SBATCH --nodes=1
#
# Расчетное время, после истечения которого задача будет принудительно
# остановлена. В данном случае --- 7 дней.
#SBATCH --time=14-00:00:00
#
# Количество потоков одного процессора (20 для intelv3-batch, 24 для
# intelv4-batch, 256 для knl-batch).
#SBATCH --ntasks-per-node=8

cd ~/Science/2018_Melatonone/DMPCMelCG/charmm-gui/gromacs/
./README
#./rundsf.sh

#g++ -c -Ialglib/cpp/src/ -O2 ./alglib/cpp/src/*.cpp
#cp ./*.o ./alglib/obj/
#rm -f *.o
#g++ ./alglib/obj/*.o ./KKEqSolitonOMP.cpp -Ialglib/cpp/src/ -fopenmp -O2
#python run.py
