#!/bin/bash
#SBATCH -o outfile_auxilary-%J
#SBATCH -e errfile_auxilary-%J
#SBATCH -n 16
#SBATCH --mem=400GB
#SBATCH -N 1
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=xxx@xxx

module purge
module load baypass/2.4

INP=/home/lazic@ad.vti.bund.de/P1.Genomic.offset/baypass/1.intermediate.files
# Aux model to get Bayes Factors, using the median omega file
#Three repetitions to account for stochasticity, final resultt was median of the output 
g_baypass -gfile $INP/baypass.gen.input -efile $INP/baypass.environment.data.txt -auxmodel -omegafile $INP/median_mat_omega.out -outprefix aux_num1 -seed 425631 -nthreads 16
g_baypass -gfile $INP/baypass.gen.input -efile $INP/baypass.environment.data.txt -auxmodel -omegafile $INP/median_mat_omega.out -outprefix aux_num2 -seed 782 -nthreads 16
g_baypass -gfile $INP/baypass.gen.input -efile $INP/baypass.environment.data.txt -auxmodel -omegafile $INP/median_mat_omega.out -outprefix aux_num3 -seed 30256 -nthreads 16
