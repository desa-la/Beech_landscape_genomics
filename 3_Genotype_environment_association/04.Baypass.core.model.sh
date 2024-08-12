#!/bin/bash
#SBATCH -o outfile_core-%J
#SBATCH -e errfile_core-%J
#SBATCH -n 24
#SBATCH --mem=400GB
#SBATCH -N 1
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=xxx@xxx

#Three repetitions to account for stochasticity, the final result was median of the output
g_baypass -gfile baypass.gen.input -outprefix core_emp_1 -seed 13156 -nthreads 24
g_baypass -gfile baypass.gen.input -outprefix core_emp_2 -seed 71754 -nthreads 24
g_baypass -gfile baypass.gen.input -outprefix core_emp_3 -seed 5990 -nthreads 24
