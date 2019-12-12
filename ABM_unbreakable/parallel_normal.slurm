#!/bin/bash
#SBATCH -J unbreakable
#SBATCH -o unbreakable.o%j
#SBATCH -e unbreakable.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 05:10:00
#SBATCH --mail-user=calebphillips@utexas.edu
#SBATCH --mail-type=all
mpiexec ./model_code -ksp_type gmres -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 10
