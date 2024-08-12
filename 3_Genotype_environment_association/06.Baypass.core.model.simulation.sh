#### This simulation analysis is used to decide on the significance threshold
#!/bin/bash
#SBATCH -o outfile_simulation-core-%J
#SBATCH -e errfile_simulation-core-%J
#SBATCH -n 16
#SBATCH --mem=10GB
#SBATCH -N 1
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=xxx@xxx


module purge
module load baypass/2.4

INP=/home/lazic@ad.vti.bund.de/P1.Genomic.offset/baypass/1.intermediate.files

g_baypass -npop 98 -gfile $INP/G.simupods2 -nthreads 24 -outprefix anapod2
