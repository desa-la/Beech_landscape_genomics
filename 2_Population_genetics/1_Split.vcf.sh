#!/bin/bash
#SBATCH -o outfile_split-vcf-%J
#SBATCH -e errfile_split-vcf-%J
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=xxx@xxx

module purge
module load vcftools/0.1.16
module load bcftools/1.18
module load htslib/1.18

VCF=/home/lazic@ad.vti.bund.de/P1.Genomic.offset/Fst/1.itermediate.files/final.filtered.bim.indep.vcf.gz
OUT=/home/lazic@ad.vti.bund.de/P1.Genomic.offset/Fst/1.itermediate.files

vcftools --gzvcf $VCF --positions $OUT/bioclim6.adaptive.qval.snps.txt --recode --stdout | bgzip -c -@ 24 > $OUT/bioclim6.adaptive.qval.snps.final.filtered.bim.indep.vcf.gz
bcftools index $OUT/bioclim6.adaptive.qval.snps.final.filtered.bim.indep.vcf.gz
vcftools --gzvcf $VCF --exclude-positions $OUT/bioclim6.adaptive.qval.snps.txt --recode --stdout | bgzip -c -@ 24 > $OUT/bioclim6.neutral.snps.final.filtered.bim.indep.vcf.gz
bcftools index $OUT/bioclim6.neutral.snps.final.filtered.bim.indep.vcf.gz
