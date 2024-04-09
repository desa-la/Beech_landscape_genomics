#!/bin/bash
#SBATCH -o outfile_vcftools_pi-%J
#SBATCH -e errfile_vcftools_pi-%J
#SBATCH -n 3
#SBATCH -N 1
#SBATCH --partition=small-jobs
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=xxx@xxx

module purge
module load vcftools/0.1.16
module load bcftools/1.18
module load samtools/1.19.2
module load htslib/1.18

VCF=/home/lazic@ad.vti.bund.de/P1.Genomic.offset/pixy/1.intermediate.files/maf.ld.6.new.plus.monomorphic.selected.samples.altdot.vcf.gz
POP=/home/lazic@ad.vti.bund.de/P1.Genomic.offset/pixy/1.intermediate.files

#Separating East West Center
bcftools view -S $POP/East.txt -O z -o East.vcf.gz $VCF
bcftools view -S $POP/West.txt -O z -o West.vcf.gz $VCF
bcftools view -S $POP/Center.txt -O z -o Center.vcf.gz $VCF

#indexing
bcftools index East.vcf.gz
bcftools index West.vcf.gz
bcftools index Center.vcf.gz

##Calculating pi
vcftools --gzvcf East.vcf.gz --window-pi 10000 --window-pi-step 10000 --out pi_East.txt
vcftools --gzvcf West.vcf.gz --window-pi 10000 --window-pi-step 10000 --out pi_West.txt
vcftools --gzvcf Center.vcf.gz --window-pi 10000 --window-pi-step 10000 --out pi_Center.txt
