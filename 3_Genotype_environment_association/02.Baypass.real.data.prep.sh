#!/bin/bash
#SBATCH -o outfile_baypass_prep-%J
#SBATCH -e errfile_baypass_prep-%J
#SBATCH -n 24
#SBATCH -N 1
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=xxx@xxx

module purge
module load vcftools/0.1.16

VCF=/home/lazic@ad.vti.bund.de/P1.Genomic.offset/Baypass/final.filtered.bim.indep.vcf.gz

list='100 101 102 103 104 108 109 10 110 111 114 115 116 117 118 11 120 124 126 127 129 12 130 132 135 136 137 138 139 13 141 142 144 145 146 14 150 161 18 20 23 24 25 26 27 28 29 2 30 36 37 38 39 40 43 44 46 48 49 51 52 53 54 55 58 59 5 61 66 67 68 69 70 71 72 73 74 75 76 77 7 80 83 84 87 88 89 8 90 91 92 93 94 95 97 98 99 9'
for sample in ${list}
do
#input data prep
vcftools --gzvcf $VCF --keep population_${sample}.corrected.txt --counts --out ${sample}
#extracting allele counts from files
awk '{print $5}' ${sample}.frq.count > baypassprep_${sample}1
awk '{print $6}' ${sample}.frq.count > baypassprep_${sample}2
#removing letter and :
sed 's/.*://' baypassprep_${sample}1 > baypassprep_${sample}1r
sed 's/.*://' baypassprep_${sample}2 > baypassprep_${sample}2r
#manually in vim remove header
sed '1d' baypassprep_${sample}1r > baypassprep_${sample}1rr
sed '1d' baypassprep_${sample}2r > baypassprep_${sample}2rr
#combining them
paste baypassprep_${sample}1rr baypassprep_${sample}2rr | column -s $'\t' -t > baypassprep_${sample}
done

#Combining all the outputs
paste baypassprep_146 baypassprep_139 baypassprep_93 baypassprep_77 baypassprep_71 baypassprep_144 baypassprep_92 baypassprep_138 baypassprep_13 baypassprep_20 baypassprep_59 baypassprep_7 baypassprep_84 baypassprep_8 baypassprep_80 baypassprep_38 baypassprep_89 baypassprep_24 baypassprep_161 baypassprep_36 baypassprep_88 baypassprep_110 baypassprep_114 baypassprep_117 baypassprep_14 baypassprep_10 baypassprep_68 baypassprep_75 baypassprep_44 baypassprep_95 baypassprep_126 baypassprep_137 baypassprep_97 baypassprep_11 baypassprep_100 baypassprep_99 baypassprep_142 baypassprep_2 baypassprep_72 baypassprep_43 baypassprep_104 baypassprep_135 baypassprep_94 baypassprep_145 baypassprep_73 baypassprep_48 baypassprep_109 baypassprep_127 baypassprep_70 baypassprep_98 baypassprep_26 baypassprep_39 baypassprep_102 baypassprep_52 baypassprep_69 baypassprep_103 baypassprep_27 baypassprep_49 baypassprep_74 baypassprep_108 baypassprep_61 baypassprep_25 baypassprep_101 baypassprep_28 baypassprep_150 baypassprep_115 baypassprep_116 baypassprep_58 baypassprep_111 baypassprep_90 baypassprep_120 baypassprep_29 baypassprep_53 baypassprep_129 baypassprep_30 baypassprep_12 baypassprep_55 baypassprep_118 baypassprep_18 baypassprep_87 baypassprep_5 baypassprep_76 baypassprep_132 baypassprep_130 baypassprep_66 baypassprep_54 baypassprep_67 baypassprep_91 baypassprep_46 baypassprep_37 baypassprep_136 baypassprep_83 baypassprep_124 baypassprep_141 baypassprep_40 baypassprep_23 baypassprep_9 baypassprep_51 | column -s $'\t' -t > baypass.gen.input
